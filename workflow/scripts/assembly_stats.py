#!/usr/bin/env python3

import sys,os,argparse,tempfile,operator,subprocess,pathlib,contextlib
from collections import defaultdict
from Bio import SeqIO


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# input: list of show-coords hits w.r.t. the same ref and qry
# assumptions: no inclusions, only overlaps allowed
def aggr_hits(lst):
    sumWeights=0
    sumIdy=0
    lastBase=0
    qryAlnBases=0
    qryBases=0
    for hit in sorted(lst,key=lambda x:int(x[2])): # process list sorted by coordinate on query
        qBeg=int(hit[2])-1
        qEnd=int(hit[3])
        weight=int(hit[4])+int(hit[5])
        sumWeights+=weight
        sumIdy+=float(hit[6])*weight
        qryAlnBases+= (qEnd-qBeg) if qBeg >= lastBase else (qEnd-lastBase)
        lastBase=qEnd # last position covered in query
        qryBases=int(hit[8])
    avgIdy = sumIdy/float(sumWeights) if sumWeights!=0 else 0.0
    covBases = (100.0*qryAlnBases)/qryBases if qryBases!=0 else 0.0
    return [avgIdy,covBases,qryAlnBases]

def fasta_n50(path,genome_size=0):
    seqLengths=[]
    for rec in SeqIO.parse(path,'fasta'):
        seqLengths.append(len(rec.seq))
    seqLengths.sort(reverse=True)

    ref_size = genome_size if genome_size > 0 else sum(seqLengths)
    n50=0
    sum_len=0
    for seq_len in seqLengths:
        sum_len+=seq_len
        if sum_len >= 0.5 * ref_size:
            n50=seq_len
            break
    return int(n50)

def parse_dnadiff_report(path):
    reportDict={}
    with open(path,'r') as report_file:
        parsed_1to1=False
        for line in report_file:
            cols=line.strip().split()
            if len(cols) > 0 and cols[0] not in reportDict:
                reportDict[cols[0]] = cols
    return reportDict

def dnadiff_dup_bases(path):
    dup_bases=0
    with open(path,'r') as diff_file:
        for line in diff_file:
            cols=line.split('\t')
            if cols[1] == 'DUP':
                dup_bases+=int(cols[4])
    return dup_bases


def main( argv = None ):

    # GET PARAMETERS
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--fasta', dest='fastaFile', required=True, help='FASTA file to be mapped against references')
    parser.add_argument('-r','--references', dest='refcsvFile', required=True, help='References CSV file: ref_id, ref_fasta_path')
    parser.add_argument('-o','--output-dir', dest='outDir', required=True, help='Output directory of results')
    parser.add_argument('-t','--threads', dest='nthreads', type=int, default=4, help='Number of threads to use with mummer')
    opt = parser.parse_args()

    # VALIDATE PARAMETERS
    validParameters = True
    if not os.path.isfile(opt.fastaFile):
        validParameters = False
        eprint(f'-f|--fasta file \"{opt.fastaFile}\" does not exist.')
    if not os.path.isfile(opt.refcsvFile):
        validParameters = False
        eprint(f'-r|--references file \"{opt.refcsvFile}\" does not exist.')
    if not validParameters:
        return 1
    
    # read reference csv file
    refSeqDict={}
    with open(opt.refcsvFile,'r') as refcsv:
        for line in refcsv:
            ref_id,ref_fasta=line.rstrip().split(',')
            refSeqDict[ref_id]=ref_fasta

    qrySeq={ qrec.id:qrec for qrec in SeqIO.parse(opt.fastaFile,'fasta') }
    qryLen={ qid:len(qrec) for qid,qrec in qrySeq.items() }

    # run MUMMer to align references to fasta, store results in a temporary directory
    qryMumDict = defaultdict(list)
    with tempfile.TemporaryDirectory() as tmpdir:
        eprint("mapping sequences to best reference")
        for ref_id in refSeqDict:
            # nucmer
            nucmer_cmd=['nucmer','--maxmatch','-t',f'{opt.nthreads}','-p',os.path.join(tmpdir,ref_id),refSeqDict[ref_id],opt.fastaFile]
            nucmer=subprocess.call(nucmer_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if nucmer != 0:
                eprint(f'error running command: {" ".join(nucmer_cmd)}')
                return 2
            # delta-filter
            with open(os.path.join(tmpdir,f'{ref_id}.qdelta'),'w') as deltaout:
                delta_cmd=['delta-filter','-q',os.path.join(tmpdir,f'{ref_id}.delta')]
                delta=subprocess.call(delta_cmd, stdout=deltaout, stderr=subprocess.DEVNULL)
                if delta != 0:
                    eprint(f'error running command: {" ".join(delta_cmd)}')
                    return 2
            # show-coords
            with open(os.path.join(tmpdir,f'{ref_id}.qcoords'),'w') as qcoords:
                showcoords_cmd=['show-coords','-qHTlc',os.path.join(tmpdir,f'{ref_id}.qdelta')]
                showcoords=subprocess.call(showcoords_cmd, stdout=qcoords, stderr=subprocess.DEVNULL)
                if showcoords != 0:
                    eprint(f'error running command: {" ".join(showcoords_cmd)}')
                    return 2
            # read mummer output and save hits to qryMumDict
            with open(os.path.join(tmpdir,f'{ref_id}.qcoords'),'r') as qcoords:
                for line in qcoords:
                    line = line.rstrip('\n')
                    if line.startswith('#') or len(line)==0:
                        continue
                    # [R_BEG] [R_END] [Q_BEG] [Q_END] [R_HITLEN] [Q_HITLEN] [%IDY] [R_LEN] [Q_LEN] [R_COV] [Q_COV] [R_NAME] [Q_NAME]
                    hit = line.split('\t') 
                    qry_id = hit[12]
                    if int(hit[2]) > int(hit[3]):
                        hit[2],hit[3] = hit[3],hit[2]
                    qryMumDict[(qry_id,ref_id)].append(hit)

    # aggregate hits concerning the same query
    qryHitList=defaultdict(list)
    for p,lst in qryMumDict.items():
        qry_id,ref_id = p
        qry_ani,qry_frac,qry_bases = aggr_hits(lst)
        if qry_frac >= 50:
            qryHitList[qry_id].append((ref_id,qry_ani,qry_frac,qry_bases))
    
    # create output directory if it does not exists
    pathlib.Path(opt.outDir).mkdir(parents=True,exist_ok=True)

    # assign best refrences to each query sequence
    with open(f'{opt.outDir}/assembly.bestref.tsv','w') as outfile:
        for qry,qry_len in qryLen.items():
            # assign queries to references w.r.t. alignment coverage/identity
            outrow=[qry]
            if len(qryHitList[qry]) > 0:
                qryHitList[qry].sort(key=lambda hit:hit[1]*hit[3],reverse=True)
                best_ref,best_ani,best_cov,best_alnbases=qryHitList[qry][0]
                outrow.append(best_ref)
                outrow.append(f'{best_ani:.2f}')
                outrow.append(f'{best_cov:.2f}')
            else:
                outrow += ['none','*','*']
            outrow.append(f'{qry_len}')
            outfile.write('\t'.join(outrow)+'\n')

    # output strain-specific files
    eprint("generating strain-specific fasta files according to best-reference")
    with contextlib.ExitStack() as stack:
        ref_list = list(refSeqDict.keys()) + ['none']
        ssOutFile = { ref_id:stack.enter_context(open(os.path.join(opt.outDir,f'assembly.{ref_id}.fa'),'w')) for ref_id in ref_list }
        for qry,sortedhits in qryHitList.items():
            if len(sortedhits) == 0:
                ssOutFile['none'].write(qrySeq[qry].format('fasta'))
                continue
            for hit in sortedhits:
                ref_id = hit[0]
                ssOutFile[ref_id].write(qrySeq[qry].format('fasta'))
                break # take the first hit only
    
    # run dnadiff and compute evaluation metrics
    eprint("running dnadiff and computing evaluation metrics")
    pathlib.Path( os.path.join(opt.outDir,'dnadiff') ).mkdir(parents=True,exist_ok=True)
    with open( os.path.join(opt.outDir,'report.tsv'), 'w') as outReport:
        outHeader = [ 
                '#ref_id', 'seq_num', 'ref_size', 'asm_size', 'N50', 'unaligned_ref', 'unaligned_asm', 'ANI', 
                'dup_ratio', 'dup_bases', 'cmp_bases', 'snps', 'inversions', 'relocations', 'transloc'
                ]
        outReport.write('\t'.join(outHeader) + '\n')
        for ref_id,ref_fasta in refSeqDict.items():
            asm_fasta = os.path.join(opt.outDir,f'assembly.{ref_id}.fa')
            if os.stat(asm_fasta).st_size == 0:
                eprint(f'warning: {asm_fasta} is empty')
                continue
            nucmer_cmd = [
                    'nucmer',
                    '--maxmatch',
                    '-t', f'{opt.nthreads}',
                    '-p', os.path.join(opt.outDir,'dnadiff',f'assembly.{ref_id}'),
                    ref_fasta,
                    asm_fasta, 
                    ]
            nucmer_proc = subprocess.call(nucmer_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if nucmer_proc != 0:
                eprint(f'error running command: {" ".join(nucmer_cmd)}')
                return 2
            dnadiff_cmd = [
                    'dnadiff',
                    '-p', os.path.join(opt.outDir,'dnadiff',f'assembly.{ref_id}'),
                    '-d', os.path.join(opt.outDir,'dnadiff',f'assembly.{ref_id}.delta'),
                    ]
            dnadiff_proc = subprocess.call(dnadiff_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if dnadiff_proc != 0:
                eprint(f'error running command: {" ".join(dnadiff_cmd)}')
                return 2
            report = parse_dnadiff_report(os.path.join(opt.outDir,'dnadiff',f'assembly.{ref_id}.report'))
            seq_num = int(report['TotalSeqs'][2])
            ref_size = int(report['TotalBases'][1])
            asm_size = int(report['TotalBases'][2])
            N50 = fasta_n50( os.path.join(opt.outDir,f'assembly.{ref_id}.fa'), ref_size )
            unaligned_ref_bases = int(report['UnalignedBases'][1].split('(')[0])
            unaligned_ref = 100.0 * unaligned_ref_bases / ref_size
            unaligned_asm_bases = int(report['UnalignedBases'][2].split('(')[0])
            unaligned_asm = 100.0 * unaligned_asm_bases / asm_size
            ANI = float(report['AvgIdentity'][1])
            aligned_ref_bases = int(report['AlignedBases'][1].split('(')[0])
            aligned_asm_bases = int(report['AlignedBases'][2].split('(')[0])
            dup_ratio = aligned_asm_bases / float(aligned_ref_bases)
            dup_bases = dnadiff_dup_bases( os.path.join(opt.outDir,'dnadiff',f'assembly.{ref_id}.qdiff') )
            cmp_bases = dnadiff_dup_bases( os.path.join(opt.outDir,'dnadiff',f'assembly.{ref_id}.rdiff') )
            snps = int(report['TotalSNPs'][1])
            inversions = int(report['Inversions'][2])
            relocations = int(report['Relocations'][2])
            transloc = int(report['Translocations'][2])
            outCols = [ 
                    f'{ref_id}', f'{seq_num}', f'{ref_size}', f'{asm_size}', f'{N50}', f'{unaligned_ref:.2f}', f'{unaligned_asm:.2f}', f'{ANI:.2f}', 
                    f'{dup_ratio:.2f}', f'{dup_bases}', f'{cmp_bases}', f'{snps}', f'{inversions}', f'{relocations}', f'{transloc}'
                    ]
            outReport.write('\t'.join(outCols) + '\n')

    return 0


if __name__ == "__main__":
    sys.exit(main())
