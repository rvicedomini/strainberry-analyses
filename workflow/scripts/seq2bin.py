#!/usr/bin/env python3

import os,sys,argparse
from Bio import SeqIO


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def fasta_ids(path):
    for seq in SeqIO.parse(path,'fasta'):
        yield seq.id

def main( argv = None ):

    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--bin-dir',  dest='binDir', required=True, help='Directory of bins')
    parser.add_argument('-o','--output',    dest='outFile', required=True, help='Output directory of results')
    parser.add_argument('-x','--extension', dest='binExt',  default='fa',  help='Extension of bins (other files in directory are ignored)')
    opt = parser.parse_args()

    if not os.path.isdir(opt.binDir):
        eprint(f'-i|--bin-dir directory "{opt.binDir}" does not exists')
        return 1

    with os.scandir(opt.binDir) as it, open(opt.outFile,'w') as out:
        for entry in it:
            if entry.name.endswith(f'.{opt.binExt}'):
                binid   = os.path.splitext(entry.name)[0]
                for seqid in fasta_ids(os.path.join(opt.binDir,entry.name)):
                    out.write(f'{seqid}\t{binid}\n')
    
    return 0


# Check if the program is not being imported
if __name__ == "__main__":
    sys.exit(main())
