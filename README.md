# Strainberry analysis workflows

+ [System requirements](#system-requirements)
+ [Before running the workflows](#before-running-the-workflows)
+ [Usage](#usage)
+ [Output files](#output-files)
+ [Cluster execution](#cluster-execution)

## System requirements

Strainberry analysis scripts make use of Snakemake and have been tested under a Linux environment.
They have been developed and tested using the following packages/tools (older versions might now work properly):

+ GNU bash (version 4.1.2(2)-release)
+ [miniconda3](https://conda.io/en/latest/miniconda.html) (Python 3.7)
+ [Snakemake](https://snakemake.readthedocs.io) (version 5.26.1)
+ Standard development packages:
    - Debian/Ubuntu: `build-essential` and `python3-dev`
    - RedHat/CentOS/Fedora: `gcc`, `gcc-c++`, `glibc-devel`, `make`, and `python3-devel` (a GCC version supporting C++11 is required)
+ [Strainberry](https://github.com/rvicedomini/strainberry) (version 1.0)

## Before running the workflows

Assembly and analysis scripts require the download of publicly available data.
Therefore it is assumed `snakemake` command to be executed from a computing node with internet access.
It is however possible to run locally the jobs performing the download of the data (*e.g.*, from a frontend node with internet access) 
and the rest of the computations on dedicated cluster nodes without internet access. See [Cluster execution](#cluster-execution) section
for more details.

## Usage

Download the workflow repository and change directory:

```
git clone https://github.com/rvicedomini/strainberry-analyses.git
```

The following configuration files are available in the `config` sub-directory in order to generate 
Strainberry main results on a single dataset at a time:

+ `mock3.yaml` - Mock3 dataset
+ `mock9_n2.yaml` - Mock9 dataset with a two-strain separation
+ `mock9_n3.yaml` - Mock9 dataset with a three-strain separation
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

## Output files


## Cluster execution


