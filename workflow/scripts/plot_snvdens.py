#!/usr/bin/env python3

import sys,os,argparse,operator
from collections import defaultdict
import matplotlib.pyplot as plt


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def main( argv = None ):

    # GET PARAMETERS
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', dest='inputFile', required=True, help='TSV file with sequence-id, SNV-density, SNV-count, sequence-length')
    parser.add_argument('-p','--prefix', dest='prefix', default="out", help='Prefix of output files')
    opt = parser.parse_args()

    # VALIDATE PARAMETERS
    validParameters = True
    if not os.path.isfile(opt.inputFile):
        validParameters = False
        eprint(f'-i|--input file \"{opt.inputFile}\" does not exist.')
    if not validParameters:
        return 1

    vatypica = [ 'SCPN01000144.1', 'SCPN01000522.1', 'SCPN01000102.1' ]
    eeligens = [ 'SCPN01000038.1', 'SCPN01001440.1', 'SCPN01000031.1' ]

    # read reference tsv file
    xlst=[]
    ylst=[]
    vatypica_xlst=[]
    vatypica_ylst=[]
    eeligens_xlst=[]
    eeligens_ylst=[]
    with open(opt.inputFile,'r') as tsvFile:
        for line in tsvFile:
            line=line.rstrip()
            if line == "":
                continue
            seqid,snvdens,_,seqlen=line.split('\t')
            if seqid not in vatypica+eeligens:
                xlst.append(float(snvdens))
                ylst.append(float(seqlen)/1e6)
            elif seqid in vatypica:
                vatypica_xlst.append(float(snvdens))
                vatypica_ylst.append(float(seqlen)/1e6)
            elif seqid in eeligens:
                eeligens_xlst.append(float(snvdens))
                eeligens_ylst.append(float(seqlen)/1e6)

    lathe, = plt.plot(xlst,ylst,'ro',label='Lathe')
    vatyp_pnts, = plt.plot(vatypica_xlst,vatypica_ylst,'bX',label='V. atypica')
    eelig_pnts, = plt.plot(eeligens_xlst,eeligens_ylst,'k^',label='E. eligens')
    #plt.xlim(right=6)
    plt.xlabel('SNV density (SNVs per 100 bp)')
    plt.ylabel('contig length (Mbp)')
    plt.legend(handles=[vatyp_pnts,eelig_pnts,lathe], labels=['V. atypica bin','E. eligens bin','Other lathe contigs'])
    #plt.axis([0, 6, 0, 20])
    plt.savefig(f'{opt.prefix}.pdf')

    return 0


# Check if the program is not being imported
if __name__ == "__main__":
    sys.exit(main())

