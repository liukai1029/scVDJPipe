from os import listdir
from os.path import isfile, join
import os
import logging
from utils import StreamToLogger
import sys
import subprocess 


#wd = "/fcrbiouatappn01/home/kliu6/projects/single_B_miseq/FASTQ_subset_test/"
wd = "/fcrbiouatappn01/home/kliu6/projects/single_B_miseq/FASTQ/"

logging.basicConfig(filename=f"{wd}log.txt",
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)
 

fl = [f for f in listdir(wd) if isfile(join(wd, f))]
fl = [f for f in fl if f[-8:] == "fastq.gz"]
fl = sorted(fl)
fp = [wd+f for f in fl]


while fp:
    input = fp.pop(0)
    mate1 = input.replace(".fastq.gz", "")
    if "_R1_" in mate1:
        mate2 = mate1.replace("_R1_","_R2_")
    if "_R2_" in mate1:
        mate2 = mate1.replace("_R2_","_R1_")
    logging.info(f"file_1: {mate1}")
    logging.info(f"file_2: {mate2}")
    
    if mate2+".fastq.gz" in fp:
        fp.pop(fp.index(mate2+".fastq.gz"))
        logging.info(f"pair found:\n{mate1}\n{mate2}\n")
        mate_id = mate1.replace("_R1_001", "")
        logging.info(f"/bin/bash -c \"/fcrbiouatappn01/home/kliu6/packages/bbmap/bbmerge.sh in1={mate1+'.fastq.gz'} in2={mate2+'.fastq.gz'} out={mate_id+'_merged.fastq'} outu={mate_id+'_merge_failed.fastq'} 2>&1 | tee -a {wd}stdout.txt\"")
        os.system(f"/bin/bash -c \"/fcrbiouatappn01/home/kliu6/packages/bbmap/bbmerge.sh in1={mate1+'.fastq.gz'} in2={mate2+'.fastq.gz'} out={mate_id+'_merged.fastq'} outu={mate_id+'_merge_failed.fastq'} 2>&1 | tee -a {wd}stdout.txt\"")
        #returned_text = subprocess.check_output(f"/bin/bash -c \"/fcrbiouatappn01/home/kliu6/packages/bbmap/bbmerge.sh in1={mate1+'.fastq.gz'} in2={mate2+'.fastq.gz'} out={mate_id+'_merged.fastq'} outu={mate_id+'_merge_failed.fastq'}\"", shell=True, universal_newlines=True) 
        #import pdb; pdb.set_trace()
        #returned_text = subprocess.check_output(f"/fcrbiouatappn01/home/kliu6/packages/bbmap/bbmerge.sh in1={mate1+'.fastq.gz'} in2={mate2+'.fastq.gz'} out={mate_id+'_merged.fastq'} outu={mate_id+'_merge_failed.fastq'}", shell=True, universal_newlines=True) 
        #import pdb; pdb.set_trace()
        #os.system(f"ls -l > stdout.log")
        #logging.info(returned_text)

