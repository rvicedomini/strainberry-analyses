#!/usr/bin/env python3

import sys,os,argparse,operator
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


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
        eprint(f'-A|--after file \"{opt.afterFile}\" does not exist.')
    if not os.path.isfile(opt.beforeFile):
        validParameters = False
        eprint(f'-B|--before file \"{opt.beforeFile}\" does not exist.')
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
            xbefore.append(float(seqcov))
            ybefore.append(float(seqlen)/1e6)

    # read reference tsv file
    xlst=[]
    ylst=[]
    asmSumCov=defaultdict(dict)
    asmSumLen=defaultdict(dict)
    with open(opt.afterFile,'r') as tsvFile:
        for line in tsvFile:
            line=line.rstrip()
            if line == "":
                continue
            seqid,seqlen,seqcov=line.split('\t')
            if not seqid.startswith('sberry|'):
                xlst.append(float(seqcov))
                ylst.append(float(seqlen)/1e6)
            else:
                seqid=seqid.split('|',1)[1] # remove "sberry|" prefix
                refname,psid,hapid,ctgidx=seqid.rsplit('_',3)
                if hapid not in asmSumCov[psid]:
                    asmSumCov[psid][hapid]=float(seqcov)*float(seqlen)/1e6
                    asmSumLen[psid][hapid]=float(seqlen)/1e6
                else:
                    asmSumCov[psid][hapid]+=float(seqcov)*float(seqlen)/1e6
                    asmSumLen[psid][hapid]+=float(seqlen)/1e6

    asmSingle_xlst=[]
    asmSingle_ylst=[]
    asmMin_xlst=[]
    asmMin_ylst=[]
    asmMax_xlst=[]
    asmMax_ylst=[]
    for psid in asmSumCov:
        if len(asmSumCov[psid])==2:
            len1=asmSumLen[psid]['h1']
            cov1=asmSumCov[psid]['h1']/len1
            len2=asmSumLen[psid]['h2']
            cov2=asmSumCov[psid]['h2']/len2
            asmMin_xlst.append( cov1 if cov1 < cov2 else cov2 )
            asmMin_ylst.append( len1 if cov1 < cov2 else len2 )
            asmMax_xlst.append( cov1 if cov1 >= cov2 else cov2 )
            asmMax_ylst.append( len1 if cov1 >= cov2 else len2 )
        elif len(asmSumCov[psid])==1:
            len1=asmSumLen[psid]['h1'] if 'h1' in asmSumLen[psid] else asmSumLen[psid]['h2']
            cov1=asmSumCov[psid]['h1']/len1 if 'h1' in asmSumLen[psid] else asmSumCov[psid]['h2']/len1
            asmSingle_xlst.append(cov1)
            asmSingle_ylst.append(len1)
    
    fig,ax = plt.subplots()
    
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    plt.plot(xbefore,ybefore,'ro',xlst,ylst,'bX',asmSingle_xlst,asmSingle_ylst,'bX',asmMin_xlst,asmMin_ylst,'bv',asmMax_xlst,asmMax_ylst,'b^',markersize=11)
    plt.xlabel('depth of coverage')
    plt.ylabel('contig length (Mbp)')

    plt.savefig(f'{opt.prefix}.pdf')

    return 0


# Check if the program is not being imported
if __name__ == "__main__":
    sys.exit(main())

