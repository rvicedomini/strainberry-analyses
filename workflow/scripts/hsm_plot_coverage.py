#!/usr/bin/env python3

import sys,os,argparse,operator
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def main( argv = None ):

    # GET PARAMETERS
    parser = argparse.ArgumentParser()
    parser.add_argument('-A','--after', dest='afterFile', required=True, help='TSV file with sequence-id, length, mean-coverage')
    parser.add_argument('-B','--before', dest='beforeFile', required=True, help='TSV file with sequence-id, length, mean-coverage')
    parser.add_argument('-p','--prefix', dest='prefix', default="out", help='Prefix of output files')
    opt = parser.parse_args()

    # VALIDATE PARAMETERS
    validParameters = True
    if not os.path.isfile(opt.afterFile):
        validParameters = False
        print(f'-A|--after file \"{opt.afterFile}\" does not exist.', file=sys.stderr)
    if not os.path.isfile(opt.beforeFile):
        validParameters = False
        print(f'-B|--before file \"{opt.beforeFile}\" does not exist.', file=sys.stderr)
    if not validParameters:
        return 1
    
    xbefore=[]
    ybefore=[]
    with open(opt.beforeFile,'r') as tsvFile:
        for line in tsvFile:
            line=line.rstrip()
            if not line:
                continue
            seqid,seqlen,seqcov=line.split('\t')
            if float(seqlen) >= 1000:
                xbefore.append(float(seqcov))
                ybefore.append(float(seqlen))

    # read reference tsv file
    xlst=[]
    ylst=[]
    with open(opt.afterFile,'r') as tsvFile:
        for line in tsvFile:
            line=line.rstrip()
            if line == "":
                continue
            seqid,seqlen,seqcov=line.split('\t')
            if float(seqlen) >= 1000:
                xlst.append(float(seqcov))
                ylst.append(float(seqlen))

    fig,ax = plt.subplots()
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    ax.set_yscale('log')
    plt.xlabel('depth of coverage')
    plt.ylabel('scaffold length (bp)')
    
    plt.plot( xbefore,ybefore,'ro', xlst,ylst,'bo', markersize=7 )
    plt.savefig(f'{opt.prefix}.pdf')

    return 0


# Check if the program is not being imported
if __name__ == "__main__":
    sys.exit(main())

