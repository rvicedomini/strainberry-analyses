#!/usr/bin/env python3

import sys,os,argparse,operator
from collections import defaultdict


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def main( argv = None ):

    # GET PARAMETERS
    parser = argparse.ArgumentParser()
    parser.add_argument('-k','--kraken', dest='inputFile', required=True, help='output of kraken2')
    parser.add_argument('-b','--bin', dest='binId', required=True, help='bin identifier')
    parser.add_argument('-p','--prefix', dest='prefix', default="out", help='prefix of output files')
    opt = parser.parse_args()

    # VALIDATE PARAMETERS
    validParameters = True
    if not os.path.isfile(opt.inputFile):
        validParameters = False
        eprint(f'-k|--kraken file \"{opt.inputFile}\" does not exist.')
    if not validParameters:
        return 1
    
    seqCount=0
    speciesLength=defaultdict(int)
    speciesSeqCount=defaultdict(int)
    with open(opt.inputFile,'r') as kraken:
        for line in kraken:
            seqid,species,seqlen=line.rstrip().split('\t')
            species='_'.join( species.split('(taxid',1)[0].split() )
            speciesLength[species]+=int(seqlen)
            speciesSeqCount[species]+=1
            seqCount+=1

    binLength=sum(speciesLength.values())
    for species in speciesLength:
        sp_len=speciesLength[species]
        sp_frac=sp_len/float(binLength)
        sp_nseq=speciesSeqCount[species]
        if sp_frac >= 0.1:
            print(f'{opt.binId}\t{seqCount}\t{binLength}\t{species}\t{sp_nseq}\t{sp_len}\t{sp_frac:.2f}')

    return 0


# Check if the program is not being imported
if __name__ == "__main__":
    sys.exit(main())

