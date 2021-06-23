import os
import pandas as pd
from snakemake.utils import min_version

min_version("5.18.0")

GLOBAL_REF_PATH = "/mnt/ssd/ssd_3/references"

# Folders
#
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

##### Config processing #####

sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")
print(sample_tab)

if config["lib_reverse_read_length"] == 0:
    read_pair_tags = [""]
else:
    read_pair_tags = ["_R1","_R2"]

wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),
    read_pair_tag = "(_R.)?"

##### Target rules #####
rule all:
    input:
        expand("mapped/{sample}.bam",sample = sample_tab.sample_name),
        expand("mapped/{sample}.bam.bai", sample = sample_tab.sample_name),

##### Modules #####

include: "rules/alignment_RNA.smk"