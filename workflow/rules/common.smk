from snakemake.utils import validate
#import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####
configfile: "config/config.yaml"

validate(config, schema="../schemas/config.schema.yaml")
sample = config["sample"]
nstrains = int(config["nstrains"])

