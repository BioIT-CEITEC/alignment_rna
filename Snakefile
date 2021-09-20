import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")

GLOBAL_REF_PATH = "/mnt/references/"

def fetch_data(file_path):
    if config["computing_type"] == "kubernetes":
        if isinstance(file_path, list) and len(file_path) == 1:
            return S3.remote(S3_BUCKET + "/" + file_path[0])
        else:
            return S3.remote(S3_BUCKET + "/" + file_path)
    else:
        if isinstance(file_path, list) and len(file_path) == 1:
            return file_path[0]
        else:
            return file_path

# setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference.json"),)
reference_dict = json.load(f)
f.close()
config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].values()][0]

##### Config processing #####
# Folders
#
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

# Samples
#
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")

if config["lib_reverse_read_length"] == 0:
    read_pair_tags = [""]
    paired = "SE"
else:
    read_pair_tags = ["_R1","_R2"]
    paired = "PE"

wildcard_constraints:
    sample = "|".join(sample_tab.sample_name),
    read_pair_tag = "(_R.)?",

##### Target rules #####
rule all:
    input:
        fetch_data(expand("mapped/{sample}.bam",sample = sample_tab.sample_name)),
        fetch_data(expand("mapped/{sample}.bam.bai", sample = sample_tab.sample_name)),
        fetch_data(expand("mapped/{sample}.bigWig",sample = sample_tab.sample_name)),
        fetch_data(expand("mapped/{sample}.bedGraph",sample = sample_tab.sample_name)),
        fetch_data("qc_reports/all_samples/alignment_RNA_multiqc/multiqc.html")

##### Modules #####

include: "rules/alignment_RNA.smk"