#!/usr/bin/env python3
import os
import random
import numpy
import sys
import argparse
import screed
from itertools import repeat, combinations
from concurrent.futures import ProcessPoolExecutor


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
        "--step",
        help="Step for random population size sampling. Default=10",
        type=int,
        default=10,
    )
    parser.add_argument(
        "--minpop",
        help="Minimum (starting) random population size. Default=20",
        type=int,
        default=20,
    )
    parser.add_argument(
        "--iterations",
        help="Number of iterations per population size. Default=100",
        type=int,
        default=100,
    )
    parser.add_argument(
        "--keepref",
        help="Retain reference sequence in analysis",
        default=False,
    )

    args = parser.parse_args(a)
    return args

def rarefaction(iterations, allseqsdict, subpopsize, cutoff):
    """
    Run rarefaction by counting sites for sampled
    population present in ≥Threshold% taxa
    """
    coresites_list=[]
    for iteration in iterations:
        print(f'Iteration: {iteration}')
        #Get random genome sample
        subsamplekeys = random.sample(list(allseqsdict.keys()), subpopsize)
        #Iterate each sequence and record gaps
        indexdict={}
        for key in subsamplekeys:
            seq = allseqsdict.get(key)
            pos =0
            for nuc in seq:
                if nuc in ("-", "N"):
                    if pos in indexdict:
                        indexdict[pos] = indexdict.get(pos) + 1
                    else:
                        indexdict[pos] = 1
                pos+=1
        seqlen = len(allseqsdict.get(subsamplekeys[0]))
        threshold_genomes = (round(subpopsize)/100)*cutoff
        noncoresites = 0
        for pos in indexdict:
            if indexdict.get(pos)<threshold_genomes:
                noncoresites+=1
        coresites_list.append(seqlen-noncoresites)
    return(coresites_list)


if __name__ == "__main__":

    args = add_args(sys.argv[1:])
    fn = args.alignment
    popsize = 0
    allseqsdict={}
    for record in screed.open(fn):
        if not args.keepref:
            if record.name == "Reference":
                pass
            else:
                allseqsdict[record.name]=record.sequence
                popsize += 1
        else:
            allseqsdict[record.name]=record.sequence
            popsize += 1
    percentcore = args.cutoff
    print("Finished reading in alignment")

    outname ='rarefaction.csv'
    output=open(outname, 'w')
    output.write('popsize')
    for j in range(0,args.iterations):
        output.write(f',iter{j}')
    output.write('\n')


    for i in range(args.minpop,popsize,args.step):
        popresults=[]
        listofsizes=[]
        nproc = args.nproc
        chunk_size = round(args.iterations/nproc)
        rangelist = list(range(0,args.iterations))
        chunks = [
            rangelist[i : i + chunk_size] for i in range(0, len(rangelist), chunk_size)
        ]
        subpopsize=i

        with ProcessPoolExecutor(nproc) as executor:
            results = executor.map(rarefaction, chunks, repeat(allseqsdict), repeat(subpopsize), repeat(args.cutoff) )
            for result_l in results:
                for result in result_l:
                    popresults.append(result)
        popresultstrs = [str(int) for int in popresults]
        outline= ','.join(popresultstrs)
        output.write(f'{subpopsize},{outline}\n')
        print(f"\n\n\nFinished popsize: {subpopsize}\n\n\n")
