import sys,argparse,operator
import pandas as pd
import matplotlib.cm as cm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.text as mtext

from matplotlib import rc
from matplotlib.legend_handler import HandlerPathCollection
from matplotlib.legend import Legend
import functools

from collections import defaultdict


def subtitle_decorator(handler):
    @functools.wraps(handler)
    def wrapper(legend, orig_handle, fontsize, handlebox):
        handle_marker = handler(legend, orig_handle, fontsize, handlebox)
        if handle_marker.get_alpha() == 0:
            handlebox.set_visible(False)
    return wrapper


def plot_clustered_stacked(dfall, binNames, preBinStats, postBinStats, str2col, labels=None, title=r"\textbf{Kraken2 classification pre/post bin separation}"):
    """Given a list of dataframes, with identical columns and index, create a clustered stacked bar plot. 
labels is a list of the names of the dataframe, used for the legend
title is a string for the title of the plot
H is the hatch used for identification of the different dataframe"""

    #rc('text', usetex=True)

    #Adds our decorator to all legend handler functions
    for handler in Legend.get_default_handler_map().values():
        handler.legend_artist = subtitle_decorator(handler.legend_artist)

    n_df = len(dfall)
    n_col = len(dfall[0].columns) 
    n_ind = len(dfall[0].index)
    
    bw=0.4
    fig,axe = plt.subplots()

    for df in dfall : # for each data frame
        axe = df.plot(kind="bar",
                      width=bw,
                      linewidth=1,
                      edgecolor="black",
                      stacked=True,
                      ax=axe,
                      legend=False,
                      grid=False,
                      colormap="Set3")  # make bar plots
    
    #color_lst=[plt.cm.Pastel1(i) for i in np.linspace(0, 1, 12)]
    #print(f'color_lst: {color_lst}')
    
    h,l = axe.get_legend_handles_labels() # get the handles we want to modify
    for i in range(0, n_df * n_col, n_col): # len(h) = n_col * n_df
        df=dfall[i//n_col]
        for j, pa in enumerate(h[i:i+n_col]):
            for bin_i, rect in enumerate(pa.patches):
                if df.iloc[bin_i][j] == 0:
                    rect.set_alpha(0)
                    continue
                c=str2col[df.index[bin_i]][df.columns[j]]
                rect.set_facecolor(plt.cm.Pastel1(c))
                rect.set_x(bin_i-(bw*n_df/2)+(bw*i/n_col))
    
    plt.suptitle(title,fontsize=32)
    #axe.set_title(title)
    xlabels=[]
    for binid in df.index:
        binCompl=postBinStats[binid]['completeness']
        binCoverage=preBinStats[binid]['coverage']
        binLabel=r'{name} ({coverage}$\times$)'.format(name=' '.join(binNames[binid][0:2]),coverage=int(binCoverage))
        bin_xlab=r'\textbf{{{name}}}'.format(name=binLabel) if binCompl>=70 else binLabel
        bin_xlab=r'\color{{red}} {name}'.format(name=bin_xlab) if binCoverage<35 else bin_xlab
        xlabels.append(bin_xlab)

    #xlabels=[ binNames[binid] for binid in df.index ]
    #axe.set_xticklabels(xlabels, rotation = 45)
    plt.xticks(np.arange(n_ind),xlabels,rotation='vertical')

    plt.ylabel('Classified bases (Mbp)',labelpad=50, fontsize=28)
    axe.yaxis.set_major_locator(plt.MultipleLocator(1000000))
    #axe.yaxis.set_minor_locator(plt.MultipleLocator(500000))
    axe.yaxis.set_major_formatter(plt.FuncFormatter(lambda x,y:f'{int(x/1000000)}'))

    #axe.grid(which='minor', color='black', alpha=0.2, axis='x')
    axe.grid(which='major', color='black', alpha=0.4, axis='y', linestyle="dotted")

    # Add invisible data to add another legend
    #n=[]
    #for i in range(n_df):
    #    n.append(axe.bar(0,0,color="gray",hatch=H*i))

    subtitles=[]
    for binid in df.index:
        #binLabel=r"\textbf{{ {spName} ({binCov}$\times$) }}".format(spName=binNames[binid],binCov=int(binStats[binid]['coverage']))
        binLabel=r"\textbf{{ {spName} }}".format(spName=' '.join(binNames[binid][0:2]))
        subtitles.append( mpatches.Patch(label=binLabel, alpha=0) )
        for species in df.columns:
            if any( dfall[i].loc[binid][species] for i in range(n_df) ):
                c=str2col[binid][species]
                speciesTokens=species.split('_')
                strainLabel='str. '+' '.join(speciesTokens[2:]) if len(speciesTokens[2:]) > 0 else 'str. representative'
                speciesLabel=' '.join(speciesTokens) if binNames[binid][0:2] != speciesTokens[0:2] else strainLabel
                subtitles.append( mpatches.Patch(facecolor=plt.cm.Pastel1(c), edgecolor='black', label=r"{}".format(speciesLabel) ) )

    lgnd=axe.legend(handles=subtitles,loc=[1.03,-0.37],ncol=1)
    axe.add_artist(lgnd)

    #red_patch = mpatches.Patch(color='red', label='The red data', alpha=0)
    #l1 = axe.legend(h[:n_col]+[red_patch], l[:n_col]+['pippo'], loc=[1.02,0])
    #if labels is not None:
    #    l2 = plt.legend(n, labels, loc=[1.01, 0.1]) 
    #axe.add_artist(l1)
    return fig,axe

def load_tsv(fname):
    bin2sp=defaultdict(lambda:defaultdict(int))
    binStats=defaultdict(lambda:defaultdict(float))
    spSet=set()
    with open(fname,'r') as preFile:
        for line in preFile:
            cols=line.rstrip().split('\t')
            binId=cols[0]
            binSize=int(cols[2])
            species=cols[3]
            spSize=int(cols[5])
            bin2sp[binId][species]=spSize
            bin2sp[binId]['others']=binSize
            binStats[binId]['completeness']=float(cols[7])
            binStats[binId]['coverage']=float(cols[10])
            spSet.add(species)
    for binId in bin2sp:
        binDict=bin2sp[binId]
        binDict['others']-=sum([binDict[sp] for sp in binDict if sp!='others'])
    return bin2sp,spSet,binStats

def main( argv = None ):

    # GET PARAMETERS
    parser = argparse.ArgumentParser()
    parser.add_argument('--pre', dest='preFile', required=True, help='TSV file of stats for the pre-separation bins')
    parser.add_argument('--sep', dest='sepFile', required=True, help='TSV file of stats for the post-separation bins')
    parser.add_argument('-p','--prefix', dest='prefix', required=True, help='output file prefix')
    opt = parser.parse_args()

    preStats,preSpSet,preBinStats=load_tsv(opt.preFile)
    sepStats,sepSpSet,sepBinStats=load_tsv(opt.sepFile)

    binList=list(preStats.keys())
    binNames=dict( (binid, max(preStats[binid].items(), key=operator.itemgetter(1))[0].split('_')) for binid in preStats.keys() )
    for binid in binNames:
        #binCoverage=preBinStats[binid]['coverage']
        binNames[binid]=[ x for i,x in enumerate(binNames[binid]) if i!=2 or x not in ['str.','str','strain'] ]

    spList=sorted(list(preSpSet|sepSpSet))
    spList.append('others')
    #spNames=list(map( lambda x:' '.join(x.split('_')),spList))

    strain2color=defaultdict(lambda:defaultdict(int))
    for bid in binList:
        blst=[ species for species in spList if preStats[bid][species] > 0 or sepStats[bid][species] > 0 ]
        for i,sp in enumerate(blst):
            strain2color[bid][sp]=i if sp!='others' else 11

    preDF = pd.DataFrame( np.array([ [ preStats[binId][species] for species in spList ] for binId in binList ]), index=binList, columns=spList )
    sepDF = pd.DataFrame( np.array([ [ sepStats[binId][species] for species in spList ] for binId in binList ]), index=binList, columns=spList )

    #plt.style.use('seaborn')
    tex_fonts = {
            'figure.figsize' : [23,13],
            'font.family' : 'sans-serif',
            'font.size' : 20,
            'legend.fontsize': 13,
            'text.usetex' : True,
            'text.latex.preamble': [
                r'\usepackage{xcolor}'
            ]
    }
    plt.rcParams.update(tex_fonts)

    fig,ax=plot_clustered_stacked([preDF,sepDF],binNames,preBinStats,sepBinStats,strain2color)

    ax.xaxis.set_tick_params(which=u'both',length=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #ax.set_frame_on(False)
    
    plt.subplots_adjust(bottom=0.27,left=0.1,right=0.8,top=0.93)
    plt.xticks(rotation=50, ha='right')
    plt.savefig(f'{opt.prefix}.barplot.svg')
    #plt.tight_layout()
    #plt.show()

# Check if the program is not being imported
if __name__ == "__main__":
    sys.exit(main())
