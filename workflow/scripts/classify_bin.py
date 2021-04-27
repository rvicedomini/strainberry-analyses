#!/usr/bin/env python3

import sys,os,argparse,operator
from collections import defaultdict


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def main( argv = None ):

    # GET PARAMETERS
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--kraken', dest='krakenFile', required=True, help='output of kraken2')
    parser.add_argument('-s','--seq2bin', dest='seq2binFile', required=True, help='sequence-bin mapping')
    parser.add_argument('-o','--output', dest='outFile', required=True, help='output file')
    opt = parser.parse_args()

    # VALIDATE PARAMETERS
    validParameters = True
    if not os.path.isfile(opt.krakenFile):
        validParameters = False
        eprint(f'-i|--kraken file \"{opt.krakenFile}\" does not exist.')
    if not os.path.isfile(opt.seq2binFile):
        validParameters = False
        eprint(f'-s|--seq2bin file \"{opt.seq2binFile}\" does not exist.')
    if not validParameters:
        return 1

    bin2seq=defaultdict(list)
    with open(opt.seq2binFile,'r') as f:
        for line in f:
            seqid,binid=line.strip().split('\t')
            bin2seq[binid].append(seqid)

    seqSpecies={}
    seqLength={}
    with open(opt.krakenFile,'r') as kraken:
        for line in kraken:
            line=line.strip()
            if not line or line.startswith('#'):
                continue
            _,seqid,species,seqlen,_=line.split('\t')
            species='_'.join( species.split('(taxid',1)[0].split() )
            seqSpecies[seqid]=species
            seqLength[seqid]=int(seqlen)

    with open(f'{opt.outFile}','w') as out:
        header=['#bin_id','bin_nseq','bin_size','strain','strain_nseq','strain_size','strain_frac']
        out.write('\t'.join(header)+'\n')
        for binid,seqlst in bin2seq.items():
            binLength=sum(seqLength[seqid] for seqid in seqlst)
            seqCount=len(bin2seq[binid])
            speciesLength=defaultdict(int)
            speciesSeqCount=defaultdict(int)
            for seqid in seqlst:
                species=seqSpecies[seqid]
                speciesLength[species]+=seqLength[seqid]
                speciesSeqCount[species]+=1
            for species in speciesLength:
                sp_nseq=speciesSeqCount[species]
                sp_len=speciesLength[species]
                sp_frac=sp_len/float(binLength)
                if sp_frac >= 0.1:
                    record=[f'{binid}',f'{seqCount}',f'{binLength}',f'{species}',f'{sp_nseq}',f'{sp_len}',f'{sp_frac:.2f}']
                    out.write('\t'.join(record)+'\n')

    return 0


if __name__ == "__main__":
    sys.exit(main())

