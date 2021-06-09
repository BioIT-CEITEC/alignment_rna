library(GenomeInfoDb)
library(data.table)

print(getwd())

### IN PROGRESS
## prepared manually with the use of GenomeInfoDb::genomeStyles
#organism_list <- c("Arabidopsis_thaliana","Caenorhabditis_elegans","Canis_familiaris","Cyanidioschyzon_merolae","Drosophila_melanogaster","Homo_sapiens","Mus_musculus","Oryza_sativa","Populus_trichocarpa","Rattus_norvegicus","Saccharomyces_cerevisiae","Zea_mays")
#names(organism_list) <- c("arabidopsis","name1","canis_familiaris","name2","drosophila_melanogaster")

run_all <- function(args){

  input_file <- args[1]
  output_file <- args[2]
  dest_DB <- args[3]
  organism <- args[4]
  if(organism %in% c("arabidopsis","Arabidopsis thaliana","Arabidopsis_thaliana") & 
     dest_DB == "UCSC" &
     ! "UCSC" %in% names(genomeStyles("Arabidopsis thaliana"))) dest_DB <- "TAIR9"
  
  a <- Sys.time()
  print(a)
  input_DT <- as.data.table(fread(paste("zgrep","-v","^#",input_file),
                                  sep="\t",
                                  colClasses=c("character"),
                                  header=FALSE))
  print(Sys.time() - a)
  
  print(paste("# getting names for",dest_DB))
  new_names <- mapSeqlevels(unique(input_DT$V1),dest_DB)
  print(Sys.time() - a)
  
  #print(new_names)
  if (length(new_names) == 0){
    print(paste("# copying data into",output_file))
    system2("cp",c(input_file,output_file))
    print(Sys.time() - a)
  } else {
    print(paste("# replacing new names in DT"))
    input_DT$V1 <- new_names[input_DT$V1]
    print(Sys.time() - a)
    print(paste("# saving data with replaced names into",output_file))
    write.table(input_DT[!is.na(V1)], output_file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    write.table(input_DT[is.na(V1)], paste0(output_file,".chr_NA"), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    print(Sys.time() - a)
  }
}


# run as Rscript
#script_dir <- dirname(sub("--file=", "", commandArgs()[grep("--file=", commandArgs())]))
args <- commandArgs(trailingOnly = T)
run_all(args)