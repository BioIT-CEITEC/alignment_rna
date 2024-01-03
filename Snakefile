import pandas as pd
from snakemake.utils import min_version

min_version("5.18.0")

configfile: "config.json"

GLOBAL_REF_PATH = config["globalResources"]
GLOBAL_TMPD_PATH = config["globalTmpdPath"]

os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

##### BioRoot utilities #####
module BR:
    snakefile: gitlab("bioroots/bioroots_utilities", path="bioroots_utilities.smk",branch="master")
    config: config

use rule * from BR as other_*

##### Config processing #####

sample_tab = BR.load_sample()

read_pair_tags = BR.set_read_pair_qc_tags()
pair_tag = BR.set_read_pair_tags()
paired = BR.set_paired_tags()

config = BR.load_organism()

wildcard_constraints:
    sample = "|".join(sample_tab.sample_name)

##### Target rules #####
rule all:
    input:
      expand("mapped/{sample}.{ext}",sample = sample_tab.sample_name, ext = ["bam","bam.bai"] if "get_cov_tracks" in config and not config["get_cov_tracks"] else ["bam","bam.bai","bigWig","bedGraph"]),
      "qc_reports/all_samples/alignment_RNA_multiqc/multiqc.html"

##### Modules #####

include: "rules/alignment_RNA.smk"
