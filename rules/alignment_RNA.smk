


rule mark_duplicates:
    input:  bam = "mapped/{sample}.not_markDups.bam",
            bai = "mapped/{sample}.not_markDups.bam.bai"
    output: bam = "mapped/{sample}.bam",
            bai = "mapped/{sample}.bam.bai"
    log: "logs/{sample}/mark_duplicates.log"
    threads:  8
    resources:  mem = 15
    params:
            mtx = "postQC/{sample}/MarkDuplicates/{sample}.markDups_metrics.txt",
            mark_duplicates= config["mark_duplicates"],
            rmDup = config["remove_duplicates"], # allow possibility for rm duplicates true
            UMI = config["UMI"],
            umi_usage= config["umi_usage"],
            keep_not_markDups_bam= config["keep_not_markDups_bam"],
    conda: "../wrappers/mark_duplicates/env.yaml"
    script: "../wrappers/mark_duplicates/script.py"


rule alignment_RNA:
    input:
        fastqs = expand("cleaned_fastq/{sample}{read_pair_tags}.fastq.gz",read_pair_tags = read_pair_tags),
        genome = expand("{ref_dir}/seq/{ref}.fa",ref_dir=reference_directory,ref=config["reference"])[0],
        fai_ucsc = expand("{ref_dir}/seq/{ref}.fa.fai.ucsc",ref_dir=reference_directory,ref=config["reference"])[0],
        gtf = expand("{ref_dir}/annot/{ref}.gtf",ref_dir=reference_directory,ref=config["reference"])[0],
        index = expand("{ref_dir}/index/STAR/SAindex",ref_dir=reference_directory,ref=config["reference"])[0]
    output:
        bam = "mapped/{sample}.not_markDups.bam",
        bai = "mapped/{sample}.not_markDups.bam.bai",
        transcriptom_bam = "mapped/{sample}.transcriptome.bam"
    log: "logs/{sample}/alignment.log"
    threads: 40
    resources:  mem = 34
    params:
        prefix = "mapped/{sample}/{sample}",
        strandness = config["strandness"],  #"true" or "false" - STAR parameters: strandedness, affects bedGraph (wiggle) files and XS tag in BAM
        num_mismatch= 999,  # Maximum number of mismatches; set this to high number (999) to disable and to use only perc_mismatch
        perc_mismatch= config["perc_mismatch"],
        max_intron= config["max_intron"],# Default used by ENCODE is 1000000; to turn this off set 1
        max_mate_dist=1000000,# Default used by ENCODE is 1000000; For "normal" fragments 1000 should be enough but for special cases, like chimeric we should increase this
        read_len=100,# Read length from the sequencing. Illumina sometimes reports N+1 http://seqanswers.com/forums/archive/index.php/t-31154.html; in case you change this value uncomment next line as well
        organism=config["organism"],
        map_perc= config["map_perc"], #pokud je to PE, nastavit na pevno 0.66
        map_score=config["map_score"], #pokud je to PE, nastavit na pevno 0.66
        pair = config["pair"], #SE or PE
        chim = "mapped/{sample}/{sample}Chimeric.out.bam",
    conda: "../wrappers/alignment_RNA/env.yaml"
    script: "../wrappers/alignment_RNA/script.py"


rule preprocess:
    input: raw = expand("raw_fastq/{sample}{read_pair_tags}.fastq.gz",read_pair_tags = read_pair_tags),
    output: cleaned = expand("cleaned_fastq/{sample}{read_pair_tags}.fastq.gz",read_pair_tags = read_pair_tags),
    log:    "sample_logs/{sample}/pre_alignment_processing.log"
    threads:    10
    resources:  mem = 10
    params: adaptors = config["adaptors"],
            r1u = "trimmed/{sample}_R1.discarded.fastq.gz",
            r2u = "trimmed/{sample}_R2.discarded.fastq.gz",
            trim_left1 = config["trim_left1"], # Applied only if trim left is true, trimming from R1 (different for classic:0, quant:10, sense:9)
            trim_right1 = config["trim_right1"], # Applied only if trim right is true, trimming from R1; you should allow this if you want to trim the last extra base and TRIM_LE is true as RD_LENGTH is not effective
            trim_left2 = config["trim_left2"], # Applied only if trim left is true, trimming from R2 (different for classic:0, quant:?, sense:7)
            trim_right2 = config["trim_right2"], # Applied only if trim right is true, trimming from R2; you should allow this if you want to trim the last extra base and TRIM_LE is true as RD_LENGTH is not effective
            phred = "-phred33",
            leading = 3,
            trailing = 3,
            crop = 250,
            minlen = config["min_length"],
            slid_w_1 = 4,
            slid_w_2 = 5,
            trim_stats = "trimmed/{sample}.PE.trim_stats.log"
    conda:  "../wrappers/preprocess/env.yaml"
    script: "../wrappers/preprocess/script.py"

# def test_func(wildcards):
#     trim_left1 = cfg.loc[cfg[SAMPLE] == wildcards.sample,"trim_left1"].min()
#     print(trim_left1)
#     return trim_left1

# rule preprocess_SE:
#     input: raw = "raw_fastq/{sample}.fastq.gz"
#     output: clean = "cleaned_fastq/{sample}.fastq.gz",
#     log:    run = "sample_logs/{sample}/preprocess_SE.log",
#             trim = "trimmed/{sample}.PE.trim_stats.log",
#     threads:    10
#     resources:  mem = 10
#     params: adaptors = lambda wildcards: cfg.loc[cfg[SAMPLE] == wildcards.sample,"adaptors"].min(),
#             trim_left1 = lambda wildcards: cfg.loc[cfg[SAMPLE] == wildcards.sample,"trim_left1"].min(), # Applied only if trim left is true, trimming from R1 (different for classic:0, quant:10, sense:9)
#             trim_right1 = lambda wildcards: cfg.loc[cfg[SAMPLE] == wildcards.sample,"trim_right1"].min(), # Applied only if trim right is true, trimming from R1; you should allow this if you want to trim the last extra base and TRIM_LE is true as RD_LENGTH is not effective
#             phred = "-phred33",
#             leading = 3,
#             trailing = 3,
#             crop = 250,
#             minlen = lambda wildcards: cfg.loc[cfg[SAMPLE] == wildcards.sample,"min_length"].min(),
#             slid_w_1 = 4,
#             slid_w_2 = 5,
#     conda:  "../wraps/fastq2bam_RNA/preprocess_SE/env.yaml"
#     script: "../wraps/fastq2bam_RNA/preprocess_SE/script.py"