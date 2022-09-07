import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")
configfile: "config.json"

GLOBAL_REF_PATH = "/mnt/references/"
GLOBAL_TMPD_PATH = "/tmp/"

# setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference.json"),)
reference_dict = json.load(f)
f.close()
config["species_name"] = [organism_name for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
config["organism"] = config["species_name"].split(" (")[0].lower().replace(" ","_")
config["species"] = config["species_name"].split(" (")[1].replace(")","")

##### Config processing #####
# Folders
#
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

# Samples
#

sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")

if not config["is_paired"]:
    read_pair_tags = ["SE"]
    pair_tag = [""]
    paired = "SE"
else:
    read_pair_tags = ["R1","R2"]
    pair_tag = ["_R1","_R2"]
    paired = "PE"

wildcard_constraints:
    sample = "|".join(sample_tab.sample_name)

##### Target rules #####
rule all:
    input:
      expand("mapped/{sample}.{ext}",sample = sample_tab.sample_name, ext = ["bam","bam.bai"] if "get_cov_tracks" in config and not config["get_cov_tracks"] else ["bam","bam.bai","bigWig","bedGraph"]),
      "qc_reports/all_samples/alignment_RNA_multiqc/multiqc.html"

##### Modules #####

include: "rules/alignment_RNA.smk"
