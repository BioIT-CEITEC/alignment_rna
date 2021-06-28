######################################
# wrapper for rule: mark_duplicates
######################################
import os
import subprocess
from snakemake.shell import shell


f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: mark_duplicates \n##\n")
f.close()

shell.executable("/bin/bash")

command = "mkdir -p "+os.path.dirname(snakemake.output.mtx)+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

#if isinstance(snakemake.input.bam, list):
#    snakemake.input.bam = snakemake.input.bam[0]
if snakemake.params.mark_duplicates == True:

    if snakemake.params.UMI[0] == "no":

        version = str(subprocess.Popen("picard MarkDuplicates --version 2>&1",shell=True,stdout=subprocess.PIPE).communicate()[0], 'utf-8')
        f = open(snakemake.log.run, 'at')
        f.write("## VERSION: Picard "+version+"\n")
        f.close()

        command = "export LD_BIND_NOW=1"
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

        command = "picard MarkDuplicates INPUT="+snakemake.input.bam+" OUTPUT="+snakemake.output.bam+" METRICS_FILE="+snakemake.output.mtx+" REMOVE_DUPLICATES="+str(nakemake.params.rmDup)+" \
            ASSUME_SORTED=true PROGRAM_RECORD_ID=null VALIDATION_STRINGENCY=LENIENT -Xmx"+str(snakemake.resources.mem)+"g 2>> "+snakemake.log.run+" "
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

        command = "samtools index "+snakemake.output.bam
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)



        command = "rm "+snakemake.input.bam
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

        command = "rm "+snakemake.input.bam + ".bai"
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

    else:

        version = str(subprocess.Popen("samtools --version 2>&1 | grep \"samtools\" ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
        f = open(snakemake.log.run, 'at')
        f.write("## VERSION: "+version+"\n")
        f.close()

        command = "umi_tools dedup -I " + snakemake.input.bam + " -S " + snakemake.output.bam + " --log " + snakemake.output.mtx \
        + " --extract-umi-method=read_id --umi-separator='_' --method=directional --edit-distance-threshold=0 \
        --spliced-is-unique --multimapping-detection-method=NH"
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)

        command = "samtools index -@" + str(snakemake.threads) + " " + snakemake.output.bam
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)
else:

    shell("mv -T " + snakemake.input.bam + " " + snakemake.output.bam)
    shell("mv -T " + snakemake.input.bai + " " + snakemake.output.bai)