#!/usr/bin/env python3
import sys
import argparse
from itertools import repeat, combinations
from concurrent.futures import ProcessPoolExecutor

import screed
import numpy as np


def add_args(a):
    """
    Parses arguments for program.
    """
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


def get_complementary_elements(mylist, idx):
    """
    Returns the sublist of `mylist` consisting of the complement of elements indexed by `idx`.
    """
    myarray = np.array(mylist)
    mask = np.full(len(mylist), False)
    mask[idx] = True
    complement = list(myarray[~mask])

    return complement


def get_pw_snps(pairs, coreseqindex):
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
                [i for i, e in enumerate(g1seqlist) if e in ("-", "N")]
                + [i for i, e in enumerate(g2seqlist) if e in ("-", "N")]
            )
        )

        if gapindex:
            g1nuc = get_complementary_elements(g1seqlist, gapindex)
            g2nuc = get_complementary_elements(g2seqlist, gapindex)
        else:
            g1nuc = g1seqlist
            g2nuc = g2seqlist

        snp_count = len(g1nuc) - sum(x == y for x, y in zip(g1nuc, g2nuc))
        sharedseq = len(g1nuc) - len(gapindex)
        snp_dist = snp_count / len(g1nuc)
        cor_snp = snp_dist * len(g1seqlist)
        results.append((f"{g1},{g2},{snp_count},{sharedseq},{snp_dist},{cor_snp}\n"))
        # print(f'{g1},{g2},{snp_count},{sharedseq},{fsnp}\n')

    return results


def get_core(filename, percentcore, popsize):
    indexdict = {}
    chunklist = list(range(0, 2 * popsize + 1, 2))

    for coord in chunklist:

        with open(filename) as filehandle:
            key = filehandle.readlines()[coord].rstrip()
            filehandle.seek(0)
            seq = list(filehandle.readlines()[coord + 1].rstrip())

        allpos = [i for i, val in enumerate(seq) if val in ("-", "N")]

        for pos in allpos:
            if pos in indexdict:
                indexdict[pos] = indexdict.get(pos) + 1
            else:
                indexdict[pos] = 1

        print(f"Sequence {key} added to core.")

    threshold = round((popsize / 100) * percentcore)
    cutoff = popsize - threshold
    thresholdgapindex = [int(key) for key in indexdict if indexdict.get(key) > cutoff]

    indexdict.clear()

    return thresholdgapindex


def remove_noncore(filename, thresholdgapindex, popsize):
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
        with open(filename) as filehandle:
            key = filehandle.readlines()[coord].rstrip()
            filehandle.seek(0)
            seq = list(filehandle.readlines()[coord + 1].rstrip())

        coreseqindex[key] = counter
        coreseq = "".join(get_complementary_elements(seq, thresholdgapindex))
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
        if not args.keepref:
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
    coreseqindex = remove_noncore(outname, gapindex, popsize)

    if args.nodists:
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
            results = executor.map(get_pw_snps, chunks, repeat(coreseqindex))

            for snps in results:
                for snp_line in snps:
                    pwoutput.write(snp_line)
