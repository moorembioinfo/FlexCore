#!/usr/bin/env python3
import sys
import argparse
from itertools import repeat, combinations
from concurrent.futures import ProcessPoolExecutor

import screed
import numpy as np


def add_args(a):
    parser = argparse.ArgumentParser(description=""" Test description """)
    parser.add_argument(
        "--alignment",
        "-a",
        help="Provide path and filename of alignment",
        required=True,
    )
    parser.add_argument(
        "--cutoff",
        "-c",
        help="Per-site percent core (integer). Default=95(%)",
        type=int,
        default=95,
    )
    parser.add_argument(
        "--nproc", "-p", help="Number of processes. Default: 1. ", type=int, default=1
    )
    parser.add_argument(
        "--nodists",
        help="Dont calculate SNP distances, only output core alignment. Default: False",
        required=False,
        action="store_true",
    )
    parser.add_argument(
        "--keepref",
        help="Retain the reference sequence in the core calculation and SNP distances. Default: False",
        required=False,
        action="store_true",
    )

    args = parser.parse_args(a)
    return args


def delete_noncore(mylist, idx):
    """
    Returns the sublist of `mylist` consisting of the complement of elements indexed by `idx`.
    """
    myarray = np.array(mylist)
    mask = np.full(len(mylist), True, dtype=bool)
    mask[idx] = False
    complement = list(myarray[mask])
    return complement


def get_pw_SNPs(pairs, popsize, coreseqindex):
    results = []
    print("Running pairwise...")

    for pair in pairs:
        g1, g2 = pair.split(",")
        g1line = coreseqindex.get(g1)
        g2line = coreseqindex.get(g2)

        with open("Coresites.fasta") as filehandle:
            g1seqlist = list(filehandle.readlines()[g1line].rstrip())
            filehandle.seek(0)
            g2seqlist = list(filehandle.readlines()[g2line].rstrip())

        # get sites missing in either seq of pair
        gapindex = list(
            set(
                [i for i, e in enumerate(g1seqlist) if e == "-"]
                + [i for i, e in enumerate(g2seqlist) if e == "-"]
            )
        )
        gapindex += list(
            set(
                [i for i, e in enumerate(g1seqlist) if e == "N"]
                + [i for i, e in enumerate(g2seqlist) if e == "N"]
            )
        )

        if len(gapindex) > 0:
            g1nuc = delete_noncore(g1seqlist, gapindex)
            g2nuc = delete_noncore(g2seqlist, gapindex)
        else:
            g1nuc = g1seqlist
            g2nuc = g2seqlist

        SNPcount = len(g1nuc) - sum(x == y for x, y in zip(g1nuc, g2nuc))
        sharedseq = len(g1nuc) - len(gapindex)
        SNPdist = SNPcount / len(g1nuc)
        corSNP = SNPdist * len(g1seqlist)
        results.append((f"{g1},{g2},{SNPcount},{sharedseq},{SNPdist},{corSNP}\n"))
        # print(f'{g1},{g2},{SNPcount},{sharedseq},{fSNP}\n')

    return results


def get_core(filename, percentcore, popsize):
    indexdict = {}
    chunklist = []

    for j in range(0, (popsize * 2) + 1, 2):
        chunklist.append(j)

    for coord in chunklist:
        seq = ""
        key = ""

        with open(filename) as filehandle:
            key = filehandle.readlines()[coord].rstrip()
            filehandle.seek(0)
            seq = list(filehandle.readlines()[coord + 1].rstrip())

        allpos = []
        allpos = [i for i, val in enumerate(seq) if val == "-" or val == "N"]

        for pos in allpos:
            if pos in indexdict:
                indexdict[pos] = indexdict.get(pos) + 1
            else:
                indexdict[pos] = 1

        print(f"Sequence {key} added to core")

    threshold = round((popsize / 100) * percentcore)
    cutoff = popsize - threshold
    thresholdgapindex = []

    for key in indexdict:
        presence = indexdict.get(key)

        if presence > cutoff:
            thresholdgapindex.append(int(key))

    indexdict.clear()

    return thresholdgapindex


def remove_noncore(filename, thresholdgapindex):
    coreseqindex = {}
    print(f"Number of non-core sites: {len(thresholdgapindex)}")
    print("Deleting non-core sites")
    outname = "Coresites.fasta"
    coreout = open(outname, "w")

    chunklist = []

    for j in range(0, (popsize * 2) + 1, 2):
        chunklist.append(j)

    counter = 1

    for coord in chunklist:
        seq = ""
        key = ""

        with open(filename) as filehandle:
            key = filehandle.readlines()[coord].rstrip()
            filehandle.seek(0)
            seq = list(filehandle.readlines()[coord + 1].rstrip())

        coreseqindex[key] = counter
        coreseq = "".join(delete_noncore(seq, thresholdgapindex))
        coreout.write(key + "\n")
        coreout.write(coreseq + "\n")
        counter += 2
        print(f"{key} core seq output")

    return coreseqindex


if __name__ == "__main__":

    args = add_args(sys.argv[1:])
    fn = args.alignment
    outname = f"{fn}.1l"
    output = open(outname, "w")
    popsize = 0

    for record in screed.open(fn):
        if args.keepref == False:
            if record.name == "Reference":
                pass
            else:
                output.write(">" + record.name + "\n")
                output.write(record.sequence + "\n")
                popsize += 1
        else:
            output.write(">" + record.name + "\n")
            output.write(record.sequence + "\n")
            popsize += 1

    percentcore = args.cutoff
    output.close()
    print("Finished alignment format conversion")
    print(f"Processing core genome for {popsize} genomes")

    # Run get-core and output core alignment
    gapindex = get_core(outname, percentcore, popsize)
    coreseqindex = remove_noncore(outname, gapindex)

    if args.nodists == True:
        print("Finished")
        exit()
    else:
        print("Running pairwise core sequence comparisons")
        pairslist = [
            ",".join(map(str, comb)) for comb in combinations(coreseqindex.keys(), 2)
        ]
        pwoutname = f"rSNP{percentcore}.csv"
        pwoutput = open(pwoutname, "w")
        pwoutput.write(f"g1,g2,SNPs,sharedseq,SNPdistance,adjustedSNPs\n")
        nproc = args.nproc
        chunk_size = int(len(pairslist) / nproc)
        chunks = []

        for i in range(0, len(pairslist), chunk_size):
            chunks.append(pairslist[i : i + chunk_size])

        with ProcessPoolExecutor(nproc) as executor:
            results = executor.map(
                get_pw_SNPs, chunks, repeat(popsize), repeat(coreseqindex)
            )

            for SNPs in results:
                for SNPline in SNPs:
                    pwoutput.write(SNPline)
