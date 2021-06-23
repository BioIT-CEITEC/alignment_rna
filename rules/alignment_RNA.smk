


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


def alignment_RNA_input(wildcards):
    preprocessed = "cleaned_fastq"
    if config["lib_reverse_read_length"] == 0:
        return os.path.join(preprocessed,"{sample}.fastq.gz")
    else:
        return [os.path.join(preprocessed,"{sample}_R1.fastq.gz"),os.path.join(preprocessed,"{sample}_R2.fastq.gz")]

rule alignment_RNA:
    input:
        fastqs = alignment_RNA_input,
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
