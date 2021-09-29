######################################
# wrapper for rule: alignment_RNA_multiqc
######################################
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: alignment_RNA_multiqc \n##\n")
f.close()

version = str(subprocess.Popen("conda list ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

multiqc_search_paths = "./qc_reports/*/MarkDuplicates/*" + " ./mapped/*" + " ./mapped/*/*_STARgenome/*" + " ./qc_reports/*/trimmomatic/*"

command = "multiqc -f -n " + snakemake.output.html + " " + multiqc_search_paths + \
              " --cl_config \"{{read_count_multiplier: 0.001, read_count_prefix: 'K', read_count_desc: 'thousands' }}\" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

