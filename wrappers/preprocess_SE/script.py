######################################
# wrapper for rule: preprocess_SE
######################################
import os
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)


f = open(log_filename, 'a+')
f.write("\n##\n## RULE: preprocess_SE \n##\n")
f.close()

version = str(subprocess.Popen("trimmomatic -version 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## VERSION: trimmomatic "+version+"\n")
f.close()

if int(snakemake.params.trim_left1) > 0 or int(snakemake.params.trim_right1) > 0:
  # the analysis needs explicit trimming
  version = str(subprocess.Popen("seqtk 2>&1 | grep \"[Vv]ersion\" | cut -f 2- -d \" \" ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
  f = open(log_filename, 'at')
  f.write("## VERSION: seqtk "+version+"\n")
  f.close()

  extra_flags_seqtk =  " -b "+str(snakemake.params.trim_left1) if int(snakemake.params.trim_left1) > 0 else ""
  extra_flags_seqtk += " -e "+str(snakemake.params.trim_right1) if int(snakemake.params.trim_right1) > 0 else ""

  command = "seqtk trimfq "+extra_flags_seqtk+" "+snakemake.input.raw+" | gzip -c > "+snakemake.input.raw.replace(".gz",".trim.gz")
  f = open(log_filename, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)

simpleClipThreshold = 10
# TODO: check for better settings (see: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf starting at page 5, or http://www.usadellab.org/cms/?page=trimmomatic)
if snakemake.params.adaptors != "null":
    adaptors = snakemake.params.adaptors.split(",")
    with open(os.path.dirname(snakemake.output.clean)+"/adapters.fa", "w") as adaptor_file:
        for i,adaptor in enumerate(adaptors):
            if adaptor.split("-")[0] == "1":
                adaptor_file.write(">adapt"+str(i)+"\n")
                adaptor_file.write(adaptor.split("-")[1]+"\n")
                if(len(adaptor.split("-")[1]) < 20):
                    simpleClipThreshold = min(len(adaptor.split("-")[1]) // 2, simpleClipThreshold)
            else:
                adaptor_file.write(">adapt"+str(i)+"\n")
                adaptor_file.write(adaptor+"\n")
                if(len(adaptor) < 20):
                    simpleClipThreshold = min(len(adaptor) // 2, simpleClipThreshold)


adaptors = "" if snakemake.params.adaptors == "null" else "ILLUMINACLIP:"+os.path.dirname(snakemake.output.clean)+"/adapters.fa:2:30:"+str(simpleClipThreshold)+":3"

trim_left = "" # HARD DEFAULT, as this is no longer in use, but might be in the future
#trim_left = "" if snakemake.params.trim_left == "" else "HEADCROP:"+snakemake.params.trim_left
input_r1 = snakemake.input.raw.replace(".gz",".trim.gz") if int(snakemake.params.trim_left1) > 0 or int(snakemake.params.trim_right1) > 0 else snakemake.input.raw


command = "trimmomatic SE -threads "+str(snakemake.threads)+" "+snakemake.params.phred+" "+input_r1+" "+snakemake.output.clean+"\
    "+adaptors+" CROP:"+str(snakemake.params.crop)+" LEADING:"+str(snakemake.params.leading)+"\
    TRAILING:"+str(snakemake.params.trailing)+" SLIDINGWINDOW:"+str(snakemake.params.slid_w_1)+":"+str(snakemake.params.slid_w_2)+" MINLEN:"+str(snakemake.params.minlen)+"\
    "+trim_left+" -Xmx"+str(snakemake.resources.mem)+"g >> "+snakemake.log.trim+" 2>&1 "

f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

if int(snakemake.params.trim_left1) > 0 or int(snakemake.params.trim_right1) > 0:
  command = "rm "+snakemake.input.raw.replace(".gz",".trim.gz")+" >> "+log_filename+" 2>&1 "
  f = open(log_filename, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)
