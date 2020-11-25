import sys,argparse,operator,gzip
import numpy as np
import matplotlib.pyplot as plt

from collections import defaultdict

def load_bed(fname):
    outList=[]
    with gzip.open(fname,'rt') as bedFile:
        for line in bedFile:
            cols=line.rstrip().split('\t')
            chrom=cols[0]
            begin=int(cols[1])
            end=int(cols[2])
            value=float(cols[3])
            outList.append( (begin,end,value) )
    return outList

def main( argv = None ):

    # GET PARAMETERS
    parser = argparse.ArgumentParser()
    parser.add_argument('-b','--bed', dest='bedFiles', nargs='+', required=True, help='gzipped BED file(s) with coverage values')
    parser.add_argument('-m','--max', dest='maxCov', type=float, default=2.5)
    parser.add_argument('-p', '--prefix', dest='prefix', required=True, help='output prefix for generated files')
    opt = parser.parse_args()

    colors=plt.get_cmap("tab10")(np.arange(len(opt.bedFiles), dtype=int))
    print(f'cmap={plt.get_cmap("tab10")}')
    print(f'range={np.arange(len(opt.bedFiles), dtype=int)}')

    tex_fonts = {
            'figure.figsize' : [25,4],
            'axes.linewidth': 3,
            'xtick.major.width': 3,
            'ytick.major.width': 3,
            'font.family' : 'sans-serif',
            'font.size' : 28,
            'legend.fontsize': 12,
    }
    plt.rcParams.update(tex_fonts)

    fig,ax = plt.subplots()

    ax.xaxis.set_major_locator(plt.MultipleLocator(500000))
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x,y:f'{x/1000000:.1f}'))
    
    ax.set_yticks(np.array([1,2,3]))
    #ax.yaxis.set_major_locator(plt.MultipleLocator(1))
    ax.set_ylim([None,opt.maxCov])
    
    ax.margins(x=0)
    
    xarr=np.array([])
    ylst=[]
    clst=[]
    for i,bed in enumerate(opt.bedFiles):
        bedList=load_bed(bed)
        xtmp=[]
        ytmp=[]
        for t in bedList:
            beg,end,value=t
            xtmp.append(beg)
            ytmp.append(min(value,opt.maxCov))
            xtmp.append(end-1)
            ytmp.append(min(value,opt.maxCov))
        xarr=np.array(xtmp)
        ylst.append(np.array(ytmp))
        clst.append(colors[i])
 
    #plt.plot(x,y,color=colors[i])
    #plt.fill_between(x,y,alpha=0.5)
    ax.plot(xarr,ylst[0],color='black') #clst[0])
    ax.plot(xarr,ylst[1],color='black') #clst[1])
    ax.fill_between(xarr, ylst[0], ylst[1], where=(ylst[0] > ylst[1]), color=clst[0], alpha=0.5, interpolate=True)
    ax.fill_between(xarr, ylst[0], ylst[1], where=(ylst[0] <= ylst[1]), color=clst[1], alpha=0.5, interpolate=True)
    ax.fill_between(xarr, np.minimum(ylst[0],ylst[1]), 0, color='grey', alpha=0.5, interpolate=True)


    # Plot axes labels and show the plot
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.38,left=0.06)
    plt.xlabel('Reference sequence (Mbp)',labelpad=30)
    plt.ylabel('Coverage',labelpad=30)
    plt.savefig(f'{opt.prefix}.svg')
    #plt.show()


# Check if the program is not being imported
if __name__ == "__main__":
    sys.exit(main())
