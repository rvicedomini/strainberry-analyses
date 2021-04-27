#!/usr/bin/env python3

import os,sys,argparse
from collections import defaultdict
from Bio import SeqIO


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def main( argv = None ):

    parser = argparse.ArgumentParser()
    parser.add_argument('--depth',  dest='depthFile', required=True, help='depth file returned by the first stage of metabat2')
    parser.add_argument('--seq2bin',  dest='seq2bin', required=True, help='tsv file containing a mapping of sequences to bins')
    parser.add_argument('--checkm',  dest='checkmFile', required=True, help='tsv table file returned by checkm')
    parser.add_argument('--class', dest='classFile', required=True, help='classification file of bins')
    parser.add_argument('--keep-best-only', dest='keepBest', action='store_true', default=False, help='Keep only best bins (completeness-5*contamination >= 70)')
    parser.add_argument('-o','--output', dest='output', required=True, help='output file')
    opt = parser.parse_args()

    validParameters=True
    if not os.path.isfile(opt.depthFile):
        validParameters=False
        eprint(f'--depth file "{opt.depthFile}" does not exist')
    if not os.path.isfile(opt.seq2bin):
        validParameters=False
        eprint(f'--seq2bin file "{opt.seq2bin}" does not exist')
    if not os.path.isfile(opt.checkmFile):
        validParameters=False
        eprint(f'--checkm file "{opt.checkmFile}" does not exist')
    if not os.path.isfile(opt.classFile):
        validParameters=False
        eprint(f'--class file "{opt.classFile}" does not exist')
    if not validParameters:
        return 1

    seqDepth={}
    seqLength={}
    with open(opt.depthFile,'r') as df:
        has_header=True
        for line in df:
            if has_header:
                has_header=False
                continue
            cols=line.strip().split('\t')
            seqLength[cols[0]]=float(cols[1])
            seqDepth[cols[0]]=float(cols[2])

    binSeqs=defaultdict(list)
    with open(opt.seq2bin,'r') as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith('#'):
                continue
            seqid,binid=line.split('\t')
            binSeqs[binid].append(seqid)

    binDepth={}
    for binid,seqlist in binSeqs.items():
        totlen=sum(seqLength[seqid] for seqid in seqlist if seqid in seqLength)
        w_sum=sum(seqDepth[seqid]*seqLength[seqid] for seqid in seqlist if seqid in seqLength)
        binDepth[binid]=w_sum/totlen

    checkmStats={}
    with open(opt.checkmFile,'r') as checkm:
        has_header=True
        for line in checkm:
            if has_header:
                has_header=False
                continue
            cols=line.split('\t')
            checkmStats[cols[0]] = list(map(float,cols[-3:]))

    binStats=[]
    with open(opt.classFile,'r') as classFile:
        for line in classFile:
            line=line.strip()
            if not line or line.startswith('#'):
                continue
            cols=line.split('\t')
            binid=cols[0]
            compl,contam,sh=checkmStats[binid]
            binqual=compl-5*contam
            bindepth=binDepth[binid]
            if not opt.keepBest or binqual >= 70:
                cols += [ f'{compl:.2f}',f'{contam:.2f}',f'{sh:.2f}',f'{bindepth:.2f}' ]
                binStats.append(cols)
    
    # sort bin statistics by descending completeness
    binStats.sort(key=lambda rec:float(rec[7]),reverse=True)
    
    with open(opt.output,'w') as out:
        header=['#bin_id','bin_nseq','bin_length','species','species_nseq','species_length','species_frac','compl','contam','s.h.','bin_depth']
        out.write('\t'.join(header)+'\n')
        for record in binStats:
            out.write('\t'.join(record)+'\n')
    
    return 0


# Check if the program is not being imported
if __name__ == "__main__":
    sys.exit(main())
