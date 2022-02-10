######################################
# wrapper for rule: alignment_SE_RNA
######################################
import os
import subprocess
import glob
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: alignment_RNA \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = "mkdir -p "+snakemake.params.prefix+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)


# if isinstance(snakemake.input.index, list):
#     star_index_dir = snakemake.input.index[0].replace("/SAindex","")
# else:
#     star_index_dir = snakemake.input.index.replace("/SAindex","")

star_index_dir = snakemake.input.index[0].replace("/SAindex","")

f = open(log_filename, 'at')

strandness = False if snakemake.params.strandness == "unstr" else True

if strandness == True:
    extra_flags_star_motif = "" # For STAR outSAMstrandField
    extra_flags_star_wig = " --outWigStrand Stranded" # For START bedGraph
    f.write("Running as Stranded experiment \n")
else:
    extra_flags_star_motif = " --outSAMstrandField intronMotif" # For STAR outSAMstrandField
    extra_flags_star_wig = " --outWigStrand Unstranded" # For START bedGraph
    f.write("Running as Unstranded experiment \n")
f.close()


if snakemake.params.paired == "SE":
    STAR_parameters = " --chimSegmentMin 30"
else:
    STAR_parameters = " --peOverlapMMp 0.1 --chimOutJunctionFormat 1 --chimSegmentMin 12 --chimJunctionOverhangMin 12"

command = "STAR --runMode alignReads --runThreadN " + str(snakemake.threads) + \
               " --genomeDir " + star_index_dir + \
               " --readFilesIn " + str(snakemake.input.fastqs)  + \
               " --readFilesCommand zcat" + \
               " --sjdbOverhang " + str(snakemake.params.read_len) + \
               " --sjdbGTFfile " + str(snakemake.input.gtf) + \
               " --outFileNamePrefix " + snakemake.params.prefix + \
               " --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1" + \
               " --outFilterMismatchNmax " + str(snakemake.params.num_mismatch) + \
               " --outFilterMismatchNoverReadLmax 1.0" + \
               " --outFilterMismatchNoverLmax "+ str(snakemake.params.perc_mismatch) + \
               " --alignIntronMin 20 --alignIntronMax " + str(snakemake.params.max_intron) + \
               " --alignMatesGapMax " + str(snakemake.params.max_mate_dist) +\
               " --outFilterMatchNmin 0 --outFilterScoreMinOverLread " + str(snakemake.params.map_score)+" --outFilterMatchNminOverLread " + str(snakemake.params.map_perc)+ \
               " --outSAMheaderHD @HD VN:1.4 SO:coordinate"+STAR_parameters+" --chimOutType Junctions SeparateSAMold" + \
               " --outSAMattrRGline ID:"+str(snakemake.wildcards.sample)+" PL:Illumina PU:"+str(snakemake.wildcards.sample)+" SM:"+str(snakemake.wildcards.sample) + \
               " --outSAMunmapped Within --outFilterType BySJout --outSAMattributes All" + \
               extra_flags_star_motif +" --quantMode GeneCounts TranscriptomeSAM --sjdbScore 1 --twopassMode Basic " + \
               " --outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate >> "+log_filename+" 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()

shell(command)

command = "mv " + snakemake.params.prefix + "Aligned.sortedByCoord.out.bam " + snakemake.output.bam + " >> "+log_filename+" 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

if hasattr(snakemake.output, 'transcriptome_bam'):
    command = "mv " + snakemake.params.prefix + "Aligned.toTranscriptome.out.bam " + snakemake.output.transcriptome_bam + " >> "+log_filename+" 2>&1 "
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

command = "rm -r " + snakemake.params.prefix + "*pass1" + " >> "+log_filename+" 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "samtools index -@ "+str(snakemake.threads)+ " "+ snakemake.output.bam + " " + snakemake.output.bai + " >> "+log_filename+" 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "STAR --runMode inputAlignmentsFromBAM" + \
                " --inputBAMfile " + snakemake.output.bam + \
                " --outWigType bedGraph" + extra_flags_star_wig + \
                " --outFileNamePrefix " + snakemake.params.prefix+" >> "+log_filename+" 2>&1 " # --outWigReferencesPrefix chr suitable for UCSC
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

if hasattr(snakemake.output, 'transcriptome_bam'):
    convert_to_ucsc = os.path.abspath(os.path.dirname(__file__))+ "/convert_chromosome_names.R"
    for bg_file in glob.glob(snakemake.params.prefix + '*.bg'):
        chr_file = bg_file.replace(".bg","") + "_chr.bedGraph"
        if os.stat(bg_file).st_size != 0:
            # We need to change chromosome names to visualize the output in UCSC Browser http://seqanswers.com/forums/archive/index.php/t-22504.html
            command = "Rscript "+convert_to_ucsc+" "+bg_file+" "+chr_file+" UCSC " + snakemake.params.organism + " >> "+log_filename+" 2>&1"
            f = open(log_filename, 'at')
            f.write("## COMMAND: "+command+"\n")
            f.close()
            shell(command)

            # Sort to proper order
            command = "LC_COLLATE=C sort -k1,1 -k2,2n -o " + bg_file + ".tmp " + chr_file + " >> "+log_filename+" 2>&1 "
            f = open(log_filename, 'at')
            f.write("## COMMAND: "+command+"\n")
            f.close()
            shell(command)

            command = "mv " + bg_file + ".tmp " + chr_file + " >> "+log_filename+" 2>&1 "
            f = open(log_filename, 'at')
            f.write("## COMMAND: "+command+"\n")
            f.close()
            shell(command)

            if os.stat(chr_file).st_size != 0:
                command = "bedGraphToBigWig "+ chr_file + " " + snakemake.input.fai_ucsc[0] + " " +  bg_file.replace(".bg","") + "_chr.bigWig >> "+log_filename+" 2>&1 "
                f = open(log_filename, 'at')
                f.write("## COMMAND: "+command+"\n")
                f.close()
                shell(command)


    # Sort transcriptome BAMs
    # Prepare for RSEM: sort transcriptome BAM to ensure the order of the reads, to make RSEM output (not pme) deterministic
    command = "samtools sort "+ snakemake.output.transcriptome_bam + " -@ " + str(snakemake.threads) + " -o "+ snakemake.output.transcriptome_bam + " >> " + log_filename + " 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: " + command + "\n")
    f.close()
    shell(command)

    command = "samtools index -@ " + str(snakemake.threads) + " " + snakemake.output.transcriptome_bam + " >> " + log_filename + " 2>&1 "
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

    # Chimeric to BAM
    command = "samtools view -@ " + str(snakemake.threads) +" -b " + snakemake.params.prefix + "Chimeric.out.sam" + " | samtools sort -@ " + str(snakemake.threads) + " -T tmp.sort " + \
    	"-o " + snakemake.params.prefix + "Chimeric.out.bam"  + " - >> " + log_filename + " 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

    command = "rm -f " + snakemake.params.prefix + "Chimeric.out.sam" + " >> " + log_filename + " 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

    command = "samtools index -@ " + str(snakemake.threads) +" "+  snakemake.params.prefix + "Chimeric.out.bam" + " >> " + log_filename + " 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
