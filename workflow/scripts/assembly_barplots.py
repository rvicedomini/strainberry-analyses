import sys,argparse,math
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


def main( argv = None ):

    # GET PARAMETERS
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', dest='inputFile', required=True, help='TSV file of assembly stats to plot')
    parser.add_argument('-p','--prefix', dest='prefix', required=True, help='output file prefix')
    opt = parser.parse_args()

    columns=['assembler','refname','nseq','refsize','asmsize','n50','unalnref','unalnasm','ani','dupratio','dupbases','cmpbases','snps','invers','reloc','transloc']
    df = pd.read_csv(opt.inputFile, sep="\t", names=columns)
    df['dupratio'] = df['dupratio'].apply(lambda x:x-1.0)

    tex_fonts = {
            'font.family' : 'sans-serif',
            'font.size' : 12,
            'legend.fontsize': 12,
            'xtick.major.pad': 10,
    }
    plt.rcParams.update(tex_fonts)

    # AVERAGE NUCLEOTIDE IDENTITY

    plt.figure()
    g = sns.catplot(x="assembler", y="ani", hue="refname", data=df, height=6, kind="bar", palette="muted")
    ani_min = math.floor(df['ani'].min()*10-1)/10
    g.ax.set(ylim=(ani_min,100))
    g.ax.grid(axis='y',linestyle='--')
    g.ax.set_axisbelow(True)
    g.ax.xaxis.set_tick_params(which=u'both',length=0)
    g.ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x,y:f'{x:.1f}'))
    
    g.set_ylabels('')
    g.set_xlabels('')
    g._legend.remove()
    handles = g._legend_data.values()
    labels = g._legend_data.keys()
    g.fig.legend(handles=handles, labels=labels, loc='lower center', ncol=3, frameon=False)
    g.fig.subplots_adjust(top=0.85, bottom=0.25, right=0.95, left=0.13)
    
    plt.suptitle(r'Average Nucleotide Identity', fontsize = 24)
    plt.savefig(f'{opt.prefix}.ani.barplot.pdf',transparent=True)

    # DUPLICATION RATIO
    
    plt.figure()
    g = sns.catplot(x="assembler", y="dupratio", hue="refname", data=df, height=6, kind="bar", palette="muted")
    dp_min = math.floor(df['dupratio'].min()*10)/10
    dp_max = math.ceil(df['dupratio'].max()*10)/10
    g.ax.set(ylim=(dp_min,dp_max))
    g.ax.grid(axis='y',linestyle='--')
    g.ax.set_axisbelow(True)
    g.ax.xaxis.set_tick_params(which=u'both',length=0)
    g.ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x,y:f'{x+1:.2f}'))

    g.ax.spines['top'].set_visible(False)
    g.ax.spines['right'].set_visible(False)
    g.ax.spines['bottom'].set_visible(False)
    g.ax.axhline(y=0,color='k')
    
    g.set_ylabels('')
    g.set_xlabels('')
    g._legend.remove()
    handles = g._legend_data.values()
    labels = g._legend_data.keys()
    g.fig.legend(handles=handles, labels=labels, loc='lower center', ncol=3, frameon=False)
    g.fig.subplots_adjust(top=0.85, bottom=0.25, right=0.95, left=0.13)
    
    plt.suptitle(r'Duplication Ratio', fontsize = 24)
    plt.savefig(f'{opt.prefix}.dupratio.barplot.pdf',transparent=True)

    return 0


if __name__ == "__main__":
    sys.exit(main())

