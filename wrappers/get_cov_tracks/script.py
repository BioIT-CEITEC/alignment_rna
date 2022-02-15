######################################
# wrapper for rule: get_cov_tracks
######################################
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: get_cov_tracks \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()


command = "samtools view "+str(snakemake.input.bam)+" | head -20 | wc -l"
with open(log_filename, 'at') as f:
	f.write("## COMMAND: " + command + "\n")
mapped_count = str(subprocess.Popen(command,shell=True,stdout=subprocess.PIPE).communicate()[0], 'utf-8')


if int(mapped_count) >= 5:
	command = "export TMPDIR="+snakemake.params.tmpd+";"+\
	          "bamCoverage"+\
	          " -b "+snakemake.input.bam+\
	          " -bs 5 "+\
	          " -o "+snakemake.output.bg+\
	          " -of bedgraph"+\
	          " -p "+str(snakemake.threads)+\
	          " >> "+log_filename+" 2>&1"
	f = open(log_filename, 'at')
	f.write("## COMMAND: "+command+"\n")
	f.close()
	shell(command)

	command = "export TMPDIR="+snakemake.params.tmpd+";"+\
	          "bamCoverage"+\
	          " -b "+snakemake.input.bam+\
	          " -bs 5 "+\
	          " -o "+snakemake.output.bw+\
	          " -of bigwig"+\
	          " -p "+str(snakemake.threads)+\
	          " >> "+log_filename+" 2>&1"
	f = open(log_filename, 'at')
	f.write("## COMMAND: "+command+"\n")
	f.close()
	shell(command)

else:
	command = "touch " + snakemake.output.bg
	f = open(log_filename, 'at')
	f.write("## COMMAND: "+command+"\n")
	f.close()
	shell(command)
	command = "touch " + snakemake.output.bw
	f = open(log_filename, 'at')
	f.write("## COMMAND: "+command+"\n")
	f.close()
	shell(command)
