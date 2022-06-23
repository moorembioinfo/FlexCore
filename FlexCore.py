#!/usr/bin/env python3
import os, sys, argparse, random
import numpy as np
import pandas as pd
from itertools import combinations
from itertools import repeat
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor


def add_args(a):
    parser = argparse.ArgumentParser(description=''' Test description ''')
    parser.add_argument('--alignment', '-a', help='Provide path and filename of alignment', required=True)
    parser.add_argument('--cutoff', '-c', help='Per-site percent core (integer). Default=95(%)', type=int, default=95)
    parser.add_argument('--nproc', '-p', help='Number of processes. Default: 1. ', type=int, default=1)
    parser.add_argument('--nodists', help='Dont calculate SNP distances, only output core alignment. Default: False', required=False, action='store_true')
    parser.add_argument('--keepref', help='Retain the reference sequence in the core calculation and SNP distances. Default: False', required=False, action='store_true')
    args = parser.parse_args(a)
    return(args)

def get_pw_SNPs(pairs, threshold, popsize, seqdict):

    results =[]

    print('Running pairwise...')

    def delete_not_pw_shared(mylist, idx):
        myarray = np.array(mylist)
        mask = np.full(len(mylist), True, dtype=bool)
        mask[idx] = False
        complement = list(myarray[mask])
        return(complement)

    for pair in pairs:
        #print(pair)
        g1 = pair.split(',')[0]
        g2 = pair.split(',')[1]
        g1seqlist = seqdict.get(g1)
        g2seqlist = seqdict.get(g2)

        gapindex = list(set([i for i, e in enumerate(g1seqlist) if e == '-'] + [i for i, e in enumerate(g2seqlist) if e == '-']))
        gapindex += list(set([i for i, e in enumerate(g1seqlist) if e == 'N'] + [i for i, e in enumerate(g2seqlist) if e == 'N']))

        if len(gapindex)>0:
            g1nuc = delete_not_pw_shared(g1seqlist, gapindex)
            g2nuc = delete_not_pw_shared(g2seqlist, gapindex)
        else:
            g1nuc = g1seqlist
            g2nuc = g2seqlist
        SNPcount = len(g1nuc) - sum(x == y for x, y in zip(g1nuc, g2nuc))

        sharedseq = len(g1nuc)-len(gapindex)

        fSNP = SNPcount/len(g1nuc)
        results.append((f'{g1},{g2},{SNPcount},{sharedseq},{fSNP}\n'))
        #print(f'{g1},{g2},{SNPcount},{sharedseq},{fSNP}\n')
    return(results)

def get_core(allseqsdict, percentcore):
    indexdict ={}
    for key in allseqsdict:
        seq = allseqsdict.get(key)
        counter =0
        allpos=[]
        allpos = ([i for i,val in enumerate(seq) if val=='-' or val=='N'])
        for pos in allpos:
            if pos in indexdict:
                indexdict[pos] = indexdict.get(pos)+1
            else:
                indexdict[pos] = 1
        print(f'Sequence {key} finished')



    popsize = len(allseqsdict.keys())
    threshold = round((popsize/100)*percentcore)
    cutoff = (popsize-threshold)
    print(cutoff)
    thresholdgapindex = []
    for key in indexdict:
        presence = indexdict.get(key)
        if presence > cutoff:
            thresholdgapindex.append(int(key))

    indexdict.clear()
    print(f'Number of non-core sites: {len(thresholdgapindex)}')
    print(sys.getsizeof(allseqsdict))


    print('Deleting non-core sites')

    def delete_noncore(mylist, idx):
        myarray = np.array(mylist)
        mask = np.full(len(mylist), True, dtype=bool)
        mask[idx] = False
        complement = list(myarray[mask])
        return(complement)

    allseqsdictcore = {}
    for key in allseqsdict:
        seq = allseqsdict.get(key)
        #allseqsdict[key] = delete_noncore(seq, thresholdgapindex)
        allseqsdictcore[key] = delete_noncore(seq, thresholdgapindex)
        print(f'{key} done (sites deleted)')

    #return(allseqsdict)
    return(allseqsdictcore)




if __name__=='__main__':

    args = add_args(sys.argv[1:])


    fn = args.alignment
    fh = open(fn)
    fc = fh.readlines()
    outname = f'{fn}.1l'
    output = open(outname, 'w')
    for line in fc:
        if line.startswith('>'):
            output.write('\n'+line)
        else:
            output.write(line.rstrip())
    percentcore = args.cutoff
    output.close()
    fc.clear()
    print('Finished alignment format conversion, to one line per sequence')

    fh1l = open(outname)
    fc1l = fh1l.readlines()
    #fc1l = fc1l[0:21]

    #Make sequence dictionary: name[seq]
    allcols = []
    allseqs = []
    allseqsdict ={}
    for line in fc1l:
        if line.startswith('>'):
            lineindex = fc1l.index(line)
            line = line.rstrip()
            if args.keepref==True:
                allcols.append(line.replace('>',''))
                nextline = fc1l[lineindex+1].rstrip()
                #allseqs.append(list(nextline))
                allseqsdict[line] = list(nextline)
            else:
                if line !='>Reference':
                    allcols.append(line.replace('>',''))
                    nextline = fc1l[lineindex+1].rstrip()
                    #allseqs.append(list(nextline))
                    allseqsdict[line] = list(nextline)
    fc1l.clear()
    fh1l.close()

    print(f'Extracting ≥{percentcore}% core sites')
    allseqsdictcore = get_core(allseqsdict, percentcore)
    allseqsdict.clear()

    if args.nodists ==True:
        print('Finished')
        exit()
    else:
        def get_all_pairs(genomeIDlist):
            allIDpairs = [",".join(map(str, comb)) for comb in combinations(genomeIDlist, 2)]
            return(allIDpairs)
        pairslist = get_all_pairs(allseqsdictcore.keys())

        print('Running pairwise')
        outname = f'rSNP{percentcore}.csv'
        output=open(outname, 'w')
        nproc = args.nproc
        chunk_size = int(len(pairslist)/nproc)
        chunks = []
        for i in range(0, len(pairslist), chunk_size):
            chunks.append(pairslist[i:i+chunk_size])

        with concurrent.futures.ProcessPoolExecutor(nproc) as executor:
            results = executor.map(get_pw_SNPs, chunks, repeat(percentcore), repeat(len(allseqsdictcore.keys())), repeat(allseqsdictcore))
            for fSNPs in results:
                for fSNPline in fSNPs:
                    output.write(fSNPline)
