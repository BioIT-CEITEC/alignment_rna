def alignment_RNA_multiqc_input(wildcards):
    input = {}
    input["bam"] = expand("mapped/{sample}.bam",sample = sample_tab.sample_name)
    # if read_pair_tags == ["SE"]:
    #     input["qc"] = expand("qc_reports/{sample}/cleaned_fastqc/SE_fastqc.html",sample = sample_tab.sample_name)
    # else:
    #     input["r1"] = expand("qc_reports/{sample}/cleaned_fastqc/R1_fastqc.html",sample=sample_tab.sample_name)
    #     input["r2"] = expand("qc_reports/{sample}/cleaned_fastqc/R2_fastqc.html",sample=sample_tab.sample_name)
    return input

rule alignment_RNA_multiqc:
    input:  unpack(alignment_RNA_multiqc_input)
    output: html = "qc_reports/all_samples/alignment_RNA_multiqc/multiqc.html"
    params: mark_duplicates= config["mark_duplicates"],
    log:    "logs/all_samples/alignment_RNA_multiqc.log"
    conda: "../wrappers/alignment_RNA_multiqc/env.yaml"
    script: "../wrappers/alignment_RNA_multiqc/script.py"

rule get_cov_tracks:
    input:  bam = "mapped/{sample}.bam",
    output: bw  = "mapped/{sample}.bigWig",
            bg  = "mapped/{sample}.bedGraph",
    log:    "logs/{sample}/get_cov_tracks.log"
    threads: 4
    params: tmpd = GLOBAL_TMPD_PATH
    conda:  "../wrappers/get_cov_tracks/env.yaml"
    script: "../wrappers/get_cov_tracks/script.py"

def mark_duplicates_input(wildcards):
    input = {}
    input["bam"] = "mapped/{sample}.not_markDups.bam"
    input["bai"] = "mapped/{sample}.not_markDups.bam.bai"
    input["transcriptome_bam"] = "mapped/transcriptome/{sample}.not_markDups.transcriptome.bam"
    #     input["transcriptome_bai"] = "mapped/transcriptome/{sample}.not_markDups.transcriptome.bam.bai"
    return input

rule mark_duplicates:
    input:  unpack(mark_duplicates_input)
    output: bam = "mapped/{sample}.bam",
            bai = "mapped/{sample}.bam.bai",
            transcriptome_bam = "mapped/transcriptome/{sample}.transcriptome.bam",
    log:    "logs/{sample}/mark_duplicates.log"
    threads: 8
    resources:  mem = 15
    params: mtx = "qc_reports/{sample}/MarkDuplicates/{sample}.markDups_metrics.txt",
            mark_duplicates = config["mark_duplicates"],
            rmDup = config["remove_duplicates"], # allow possibility for rm duplicates true
            UMI = config["UMI"],
            umi_usage = config["umi_usage"],
            keep_not_markDups_bam = config["keep_not_markDups_bam"],
    conda: "../wrappers/mark_duplicates/env.yaml"
    script: "../wrappers/mark_duplicates/script.py"

def alignment_RNA_input(wildcards):
    processed = "processed_fastq"
    if read_pair_tags == ["SE"]:
        return os.path.join(processed,"{sample}.fastq.gz")
    else:
        return [os.path.join(processed,"{sample}_R1.fastq.gz"),os.path.join(processed,"{sample}_R2.fastq.gz")]

rule alignment_RNA:
    input:  fastqs = alignment_RNA_input,
            genome = config["organism_fasta"], # defined in utilities
            fai_ucsc = config["organism_ucsc"], # defined in utilities
            gtf = config["organism_gtf"], # defined in utilities
            index = config["organism_star"] # defined in utilities
    output: bam = "mapped/{sample}.not_markDups.bam",
            bai = "mapped/{sample}.not_markDups.bam.bai",
            transcriptome_bam = "mapped/transcriptome/{sample}.not_markDups.transcriptome.bam",
            # transcriptome_bai = "mapped/transcriptome/{sample}.not_markDups.transcriptome.bam.bai",
    log:    "logs/{sample}/alignment.log"
    threads: 40
    resources:  mem = 34
    params: prefix = "mapped/{sample}/{sample}",
            strandness = config["strandness"],
            num_mismatch= 999,  # Maximum number of mismatches; set this to high number (999) to disable and to use only perc_mismatch
            perc_mismatch= config["perc_mismatch"],
            max_intron= config["max_intron"],# Default used by ENCODE is 1000000; to turn this off set 1
            max_mate_dist=1000000,# Default used by ENCODE is 1000000; For "normal" fragments 1000 should be enough but for special cases, like chimeric we should increase this
            read_len=100,# Read length from the sequencing. Illumina sometimes reports N+1 http://seqanswers.com/forums/archive/index.php/t-31154.html; in case you change this value uncomment next line as well
            organism=config["organism"],
            map_perc= config["map_perc"],
            map_score=config["map_score"],
            paired = paired,
            tmpd = GLOBAL_TMPD_PATH,
    conda: "../wrappers/alignment_RNA/env.yaml"
    script: "../wrappers/alignment_RNA/script.py"
