

rule alignment_RNA_multiqc:
    input:  bam = expand("mapped/{sample}.bam",sample = sample_tab.sample_name),
    output: html = "qc_reports/all_samples/alignment_RNA_multiqc/multiqc.html"
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
    # if config["RSEM"]:
    #     input["transcriptome_bam"] = "mapped/transcriptome/{sample}.not_markDups.transcriptome.bam"
    #     input["transcriptome_bai"] = "mapped/transcriptome/{sample}.not_markDups.transcriptome.bam.bai"
    return input


rule mark_duplicates:
    input:  unpack(mark_duplicates_input)
    output: bam = "mapped/{sample}.bam",
            bai = "mapped/{sample}.bam.bai",
    log:    "logs/{sample}/mark_duplicates.log"
    threads: 8
    resources:  mem = 15
    params: mtx = "qc_reports/{sample}/MarkDuplicates/{sample}.markDups_metrics.txt",
            mark_duplicates= config["mark_duplicates"],
            rmDup = config["remove_duplicates"], # allow possibility for rm duplicates true
            UMI = config["UMI"],
            umi_usage= config["umi_usage"],
            keep_not_markDups_bam= config["keep_not_markDups_bam"],
            # RSEM = config["RSEM"],
    conda: "../wrappers/mark_duplicates/env.yaml"
    script: "../wrappers/mark_duplicates/script.py"


def alignment_RNA_input(wildcards):
    preprocessed = "cleaned_fastq"
    if read_pair_tags == [""]:
        return os.path.join(preprocessed,"{sample}.fastq.gz")
    else:
        return [os.path.join(preprocessed,"{sample}_R1.fastq.gz"),os.path.join(preprocessed,"{sample}_R2.fastq.gz")]

rule alignment_RNA:
    input:  fastqs = alignment_RNA_input,
            genome = expand("{ref_dir}/seq/{ref}.fa",ref_dir=reference_directory,ref=config["reference"]),
            fai_ucsc = expand("{ref_dir}/seq/{ref}.fa.fai.ucsc",ref_dir=reference_directory,ref=config["reference"]),
            gtf = expand("{ref_dir}/annot/{ref}.gtf",ref_dir=reference_directory,ref=config["reference"]),
            index = expand("{ref_dir}/index/STAR/SAindex",ref_dir=reference_directory,ref=config["reference"])
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

def cleaned_fastq_qc_input(wildcards):
    if wildcards.read_pair_tag == "SE":
        input_fastq_read_pair_tag = ""
    else:
        input_fastq_read_pair_tag = "_" + wildcards.read_pair_tag
    return f'cleaned_fastq/{wildcards.sample}{input_fastq_read_pair_tag}.fastq.gz'

rule cleaned_fastq_qc:
    input:  cleaned_fastq = cleaned_fastq_qc_input
    output: html = "qc_reports/{sample}/cleaned_fastqc/{read_pair_tag}_fastqc.html"
    log:    "logs/{sample}/cleaned_fastqc_{read_pair_tag}.log"
    params: extra = "--noextract --format fastq --nogroup",
    threads:  1
    conda:  "../wrappers/cleaned_fastq_qc/env.yaml"
    script: "../wrappers/cleaned_fastq_qc/script.py"

rule preprocess:
    input:  raw = expand("raw_fastq/{{sample}}{read_pair_tags}.fastq.gz",read_pair_tags = read_pair_tags),
    output: cleaned = expand("cleaned_fastq/{{sample}}{read_pair_tags}.fastq.gz",read_pair_tags = read_pair_tags),
    log:    "logs/{sample}/pre_alignment_processing.log"
    threads: 10
    resources:  mem = 10
    params: adaptors = config["trim_adapters"],
            r1u = "cleaned_fastq/trimmed/{sample}_R1.discarded.fastq.gz",
            r2u = "cleaned_fastq/trimmed/{sample}_R2.discarded.fastq.gz",
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
            trim_stats = "qc_reports/{sample}/trimmomatic/trim_stats.log"
    conda:  "../wrappers/preprocess/env.yaml"
    script: "../wrappers/preprocess/script.py"

