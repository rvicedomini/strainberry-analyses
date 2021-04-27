#!/usr/bin/env python3

import os,sys,argparse
from collections import defaultdict


def main( argv = None ):

    parser = argparse.ArgumentParser()
    parser.add_argument('--bins-dir',  dest='binsDir', required=True)
    parser.add_argument('--depth-ext',  dest='depthExt', default='depth')
    parser.add_argument('--class-ext',  dest='classExt', default='kraken2')
    parser.add_argument('--checkm',  dest='checkmFile', required=True, help='tsv table file returned by checkm')
    parser.add_argument('-o','--output', dest='output', required=True, help='output file')
    opt = parser.parse_args()

    validParameters=True
    if not os.path.isdir(opt.binsDir):
        validParameters=False
        print(f'Bins directory "{opt.binsDir}" does not exist',file=sys.stderr)
    if not validParameters:
        return 1

    # load checkm stats of bins
    checkmStats={}
    with open(opt.checkmFile,'r') as checkm:
        has_header=True
        for line in checkm:
            if has_header:
                has_header=False
                continue
            cols=line.split('\t')
            checkmStats[cols[0]] = list(map(float,cols[-3:]))

    # load kraken2 classifications of bins
    seqDepth=defaultdict(dict)
    seqSpecies=defaultdict(dict)
    seqLength=defaultdict(dict)
    with os.scandir(opt.binsDir) as it:
        for entry in it:
            if entry.is_file() and entry.name.endswith(f'.{opt.classExt}'):
                binid = os.path.splitext(entry.name)[0]
                with open(os.path.join(opt.binsDir,entry.name),'r') as kraken:
                    for line in kraken:
                        line=line.strip()
                        if not line or line.startswith('#'):
                            continue
                        _,seqid,species,seqlen,_=line.split('\t')
                        species='_'.join( species.split('(taxid',1)[0].split() )
                        seqSpecies[binid][seqid]=species
                        seqLength[binid][seqid]=int(seqlen)
                        seqDepth[binid][seqid]=0.0
    

    # load bin depth of coverage
    with os.scandir(opt.binsDir) as it:
        for entry in it:
            if entry.is_file() and entry.name.endswith(f'.{opt.depthExt}'):
                binid = os.path.splitext(entry.name)[0]
                with open(os.path.join(opt.binsDir,entry.name),'r') as df:
                    for line in df:
                        cols=line.strip().split('\t')
                        if cols[0]=='contigName':
                            continue
                        seqDepth[binid][cols[0]]=float(cols[2])

    # compute bin average depth of coverage
    binDepth={}
    for binid in seqDepth:
        totlen=sum(seqLength[binid][seqid] for seqid in seqLength[binid])
        w_sum=sum(seqDepth[binid][seqid]*seqLength[binid][seqid] for seqid in seqLength[binid])
        binDepth[binid]=w_sum/totlen if totlen > 0 else 0.0

    binStats=[]
    #header=['#bin_id','bin_nseq','bin_size','strain','strain_nseq','strain_size','strain_frac']
    for binid in binDepth:
        binLength=sum(seqLength[binid][seqid] for seqid in seqLength[binid])
        seqCount=len(seqLength[binid])
        speciesLength=defaultdict(int)
        speciesSeqCount=defaultdict(int)
        for seqid in seqLength[binid]:
            species=seqSpecies[binid][seqid]
            speciesLength[species]+=seqLength[binid][seqid]
            speciesSeqCount[species]+=1
        for species in speciesLength:
            sp_nseq=speciesSeqCount[species]
            sp_len=speciesLength[species]
            sp_frac=sp_len/float(binLength)
            if sp_frac >= 0.1:
                compl,contam,sh=checkmStats[binid]
                bindepth=binDepth[binid]
                record=[f'{binid}',f'{seqCount}',f'{binLength}',f'{species}',f'{sp_nseq}',f'{sp_len}',f'{sp_frac:.2f}',f'{compl:.2f}',f'{contam:.2f}',f'{sh:.2f}',f'{bindepth:.2f}']
                binStats.append(record)

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
