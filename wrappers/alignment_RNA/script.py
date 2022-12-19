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

command = "(time STAR" + \
           " --runMode alignReads --runThreadN " + str(snakemake.threads) + \
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
           " --outFilterMatchNmin 0 --outFilterScoreMinOverLread " + str(snakemake.params.map_score)+\
           " --outFilterMatchNminOverLread " + str(snakemake.params.map_perc)+ \
           " --outSAMheaderHD @HD VN:1.4 SO:coordinate"+\
           STAR_parameters+" --chimOutType Junctions SeparateSAMold" + \
           " --outSAMattrRGline ID:"+str(snakemake.wildcards.sample)+" PL:Illumina PU:"+str(snakemake.wildcards.sample)+" SM:"+str(snakemake.wildcards.sample) + \
           " --outSAMunmapped Within --outFilterType BySJout --outSAMattributes All" + \
           extra_flags_star_motif +" --quantMode GeneCounts TranscriptomeSAM --sjdbScore 1 --twopassMode Basic " + \
           " --outMultimapperOrder Random --outSAMtype BAM SortedByCoordinate" + \
           " --outTmpDir " + snakemake.params.tmpd + "/" + snakemake.wildcards.sample + \
           ") >> "+log_filename+" 2>&1 "
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

command = "(time samtools index -@ "+str(snakemake.threads)+ " "+ snakemake.output.bam + " " + snakemake.output.bai + ") >> "+log_filename+" 2>&1 "
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)


command = "(time STAR" + \
          " --runMode inputAlignmentsFromBAM" + \
          " --inputBAMfile " + snakemake.output.bam + \
          " --outWigType bedGraph" + extra_flags_star_wig + \
          " --outFileNamePrefix " + snakemake.params.prefix+\
          " --outTmpDir " + snakemake.params.tmpd + "/" + snakemake.wildcards.sample + \
          ") >> "+log_filename+" 2>&1 " # --outWigReferencesPrefix chr suitable for UCSC
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
            command = "(time Rscript "+convert_to_ucsc+" "+bg_file+" "+chr_file+" UCSC " + snakemake.params.organism + ") >> "+log_filename+" 2>&1"
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
                command = "(time bedGraphToBigWig "+ chr_file + " " + snakemake.input.fai_ucsc[0] +\
                          " " +  bg_file.replace(".bg","") + "_chr.bigWig"+\
                          ") >> "+log_filename+" 2>&1 "
                f = open(log_filename, 'at')
                f.write("## COMMAND: "+command+"\n")
                f.close()
                shell(command)


    # THIS IS NOT WORKING, multimapped reads are not necessarily grouped together according to rsem-sam-validator tool
    # # Sort transcriptome BAMs
    # # Prepare for RSEM: sort transcriptome BAM to ensure the order of the reads, to make RSEM output (not pme) deterministic
    # if snakemake.params.paired == "SE":
    #     command = "(time cat <(samtools view -H "+snakemake.output.transcriptom_bam+")"+\
    #               " <(samtools view -@ "+str(snakemake.threads)+\
    #               " "+snakemake.output.transcriptom_bam+\
    #               " | sort -S "+str(snakemake.resources.mem)+"G -T "+snakemake.params.tmpd+")"+\
    #               " | samtools view -@ "+str(snakemake.threads)+" -b -"+\
    #               ") > "+snakemake.output.transcriptom_bam+".tmp 2>> "+log_filename
    # else:
    #     command = "(time cat <(samtools view -H "+snakemake.output.transcriptom_bam+")"+\
    #               " <(samtools view -@ "+str(snakemake.threads)+\
    #               " "+snakemake.output.transcriptom_bam+\
    #               " | awk '{{line=$0; getline; printf \"%s %s\\n\",line,$0}}'"+\
    #               " | sort -S "+str(snakemake.resources.mem)+"G -T "+snakemake.params.tmpd+\
    #               " | tr ' ' '\\n')"+\
    #               " | samtools view -@ "+str(snakemake.threads)+" -b -"+\
    #               ") > "+snakemake.output.transcriptom_bam+".tmp 2>> "+log_filename
    # f = open(log_filename, 'at')
    # f.write("## COMMAND: "+command+"\n")
    # f.close()
    # shell(command)
    #     
    # command = "mv "+snakemake.output.transcriptom_bam+".tmp " + snakemake.output.transcriptom_bam+" >> "+log_filename+" 2>&1"
    # f = open(log_filename, 'at')
    # f.write("## COMMAND: "+command+"\n")
    # f.close()
    # shell(command)

    # Chimeric to BAM
    command = "(time samtools view -@ " + str(snakemake.threads) +\
              " -b " + snakemake.params.prefix + "Chimeric.out.sam" +\
              " | samtools sort -@ " + str(snakemake.threads) + " -T "+snakemake.params.tmpd+\
              " -o " + snakemake.params.prefix + "Chimeric.out.bam -" +\
              ") >> " + log_filename + " 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

    command = "rm -f " + snakemake.params.prefix + "Chimeric.out.sam" + " >> " + log_filename + " 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

    command = "(time samtools index -@ " + str(snakemake.threads) +\
              " "+  snakemake.params.prefix + "Chimeric.out.bam" +\
              ") >> " + log_filename + " 2>&1"
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
