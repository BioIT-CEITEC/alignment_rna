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

if snakemake.params.mark_duplicates == True:
    extra = "./qc_reports/*/MarkDuplicates/*"
else:
    extra = ""
# multiqc_search_paths =  extra + " ./qc_reports/*/cleaned_fastqc/*" + " ./mapped/*" + " ./mapped/*/*_STARgenome/*" + " ./qc_reports/*/trimmomatic/*"
multiqc_search_paths =  extra + " ./mapped/*" + " ./mapped/*/*_STARgenome/*"

command = "multiqc -f -n " + snakemake.output.html + " " + multiqc_search_paths + \
              " --cl-config \"{{read_count_multiplier: 0.001, read_count_prefix: 'K', read_count_desc: 'thousands' }}\" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

