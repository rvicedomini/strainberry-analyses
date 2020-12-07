#!/usr/bin/env python3

# splits a strainberry contig assembly into bins according to the binning 
# of its input strain-oblivious assembly

import os,sys,argparse
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def insert_newlines(inString, every=80):
    return inString if every <= 0 else '\n'.join(inString[i:i+every] for i in range(0, len(inString), every))


def main( argv = None ):

    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--strainberry', dest='sberryFile', required=True, help='Strainberry contig assembly')
    parser.add_argument('-s','--seq2bin', dest='seq2binFile', required=True, help='Mapping of strain-oblivious sequences to bins')
    parser.add_argument('-b','--bin-dir', dest='binDir', required=True, help='Output directory of bins')
    parser.add_argument('-f','--force', dest='force', action='store_true', default=False, help='Overwrite output files if bin directory exists')
    parser.add_argument('-x','--extension', dest='binExt',  default='fa',  help='Extension of bins (other files in directory are ignored)')
    opt = parser.parse_args()

    validParameters = True
    if not os.path.isfile(opt.sberryFile):
        validParameters = False
        eprint(f'-i|--strainberry file "{opt.sberryFile}" does not exists')
    if not os.path.isfile(opt.seq2binFile):
        validParameters = False
        eprint(f'-s|--seq2bin file directory "{opt.seq2binFile}" does not exists')
    if not validParameters:
        return 1

    if os.path.isdir(opt.binDir) and opt.force:
        eprint(f'warning: existing files in bin directory will be overwritten')
    if os.path.isdir(opt.binDir) and not opt.force:
        eprint(f'error: bin directory "{opt.binDir}" already exists')
        return 2

    # it assumes one sequence belong to only one bin
    ctg2bin={}
    with open(opt.seq2binFile,'r') as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith('#'):
                continue
            ctgid,binid=line.split('\t')
            ctg2bin[ctgid]=binid

    sberrySeq={}
    bin2sberry=defaultdict(list)
    for sobj in SeqIO.parse(opt.sberryFile,'fasta'):
        seqid=sobj.id
        sberrySeq[seqid]=str(sobj.seq)
        ctgid = seqid.split('|',1)[1].rsplit('_',3)[0] if seqid.startswith('sberry|') else seqid.rsplit('_',1)[0]
        binid=ctg2bin[ctgid]
        bin2sberry[binid].append(seqid)
    
    Path(opt.binDir).mkdir(parents=True,exist_ok=True)
    for binid,seqlst in bin2sberry.items():
        binpath=os.path.join(opt.binDir,f'{binid}.{opt.binExt}')
        with open(binpath,'w') as out:
            for seqid in seqlst:
                out.write(f'>{seqid}\n')
                out.write(f'{insert_newlines(sberrySeq[seqid])}\n')

    return 0


if __name__ == "__main__":
    sys.exit(main())

