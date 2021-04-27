# Strainberry analysis workflows

+ [System requirements](#system-requirements)
+ [Before running the workflows](#before-running-the-workflows)
+ [Usage](#usage)
+ [Output files](#output-files)
+ [Execution in a computing environment](#execution-in-a-computing-environment)

## System requirements

+ GNU bash (version 4+ recommended)

Strainberry analysis scripts make use of Snakemake and have been tested under a Linux environment.
They have been developed and tested using the following packages/tools (older versions might not work as intended):

+ [miniconda3](https://conda.io/en/latest/miniconda.html) (Python 3.8)
+ [Snakemake](https://snakemake.readthedocs.io) (version 6.0.5)
+ [Strainberry](https://github.com/rvicedomini/strainberry) (version 1.1)

## Before running the workflows

Assembly and analysis scripts require the download of publicly available data that will be stored in the
`resources` sub-directory. Therefore, it is assumed `snakemake` command to be executed from a computing node with internet access.

In a computing environment (*e.g.*, SLURM, SGE), it is however possible to run locally the jobs performing the download of the data 
from a frontend node (with internet access) and the rest of the computations on dedicated cluster nodes 
(without internet access). See "[Execution in a computing environment](#execution-in-a-computing-environment)" section for more details.

#### Note on the Human Stool Microbiome dataset

Running the analysis workflow related to the Human Stool Microbiome dataset will make use of the
[Maxikraken2 database](https://lomanlab.github.io/mockcommunity/mc_databases.html) (March 2019 version). 
By default it is automatically downloaded in the `resources/maxikraken2_1903_140GB` directory. 
It is however possible to specify a custom location (in case it is already available in the system)
and use that path in the analysis workflow.
In order to achieve that, it is necessary to modify the last two lines of the configuration file `config/hsm.yaml` as follows:

```
download_kraken2_db: False
kraken2_db: '/absolute/path/to/maxikraken2_1903_140GB'
```

## Usage

Clone the repository:

```
git clone https://github.com/rvicedomini/strainberry-analyses.git
```

The following configuration files are available in the `config` sub-directory in order to generate 
Strainberry main results on a single dataset at a time:

+ `mock3.yaml` - Mock3 dataset
+ `mock9.yaml` - Mock9 dataset
+ `nwc2_pacbio.yaml` - NWC2 dataset (PacBio)
+ `nwc2_ont.yaml` - NWC2 dataset (Nanopore)
+ `hsm.yaml` - Human Stool Microbiome dataset

Three Snakemake workflows are provided in the `workflow` sub-directory in order to run the specific analyses
carried out on the different datasets:

+ `Snakefile-Mock.smk` - to be used with `mock*.yaml` configuration files
+ `Snakefile-NWC2.smk` - to be used with `nwc2_*.yaml` configuration files
+ `Snakefile-HSM.smk` - to be used with `hsm.yaml` configuration file

Finally, it is assumed that `snakemake` command is run from the base directory 
of the repository. A usage example for running the analysis workflow on the 
Mock3 dataset is shown below.

#### Example - Running the analysis workflow on Mock3

```
cd strainberry-analyses
snakemake --use-conda --cores 12 --snakefile workflow/Snakefile-Mock.smk --configfile config/mock3.yaml
```
will run the main analysis workflow on Mock3, using at most 12 threads.

## Results

Assemblies, evaluation files, and graphs are generated in the `results` sub-directory in paths that
depend on the dataset processed (see [Usage](#usage)).

+ Assemblies are saved in the `results/[dataset]/assemblies` directory
+ Evaluation tables and graphs are generated in the `results/[dataset]/evaluation` directory

where `[dataset]` is either `mock3`, `mock9`, `nwc2_pacbio`, `nwc2_ont`, or `hsm`.

## Execution in a computing environment

The documentation for running snakemake in a computing environment can be found at the following links:

+ [Cluster execution](https://snakemake.readthedocs.io/en/stable/executing/cluster.html)
+ [Snakemake Profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)

