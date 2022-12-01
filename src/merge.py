from os import listdir
from os.path import isfile, join
import os
import logging
from utils import StreamToLogger
import sys

wd = "/fcrbiouatappn01/home/kliu6/projects/single_B_miseq/FASTQ_subset_test/"

logging.basicConfig(filename=f"{wd}log.txt",
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)
 

 

#logging.basicConfig(
#   level=logging.DEBUG,
#   format='%(asctime)s:%(levelname)s:%(name)s:%(message)s',
#   filename="out.log",
#   filemode='a'
#)
#
#stdout_logger = logging.getLogger('STDOUT')
#sl = StreamToLogger(stdout_logger, logging.INFO)
#sys.stdout = sl
#
#stderr_logger = logging.getLogger('STDERR')
#sl = StreamToLogger(stderr_logger, logging.ERROR)
#sys.stderr = sl
#
#




fl = [f for f in listdir(wd) if isfile(join(wd, f))]
fl = [f for f in fl if f[-8:] == "fastq.gz"]
fl = sorted(fl)
fp = [wd+f for f in fl]

#import pdb; pdb.set_trace()

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
        logging.info(f"/bin/bash -c \"/fcrbiouatappn01/home/kliu6/packages/bbmap/bbmerge.sh in1={mate1+'.fastq.gz'} in2={mate2+'.fastq.gz'} out={mate_id+'_merged.fastq'} outu={mate_id+'_merge_failed.fastq'} > {wd}stdout.txt\"")
        os.system(f"/bin/bash -c \"/fcrbiouatappn01/home/kliu6/packages/bbmap/bbmerge.sh in1={mate1+'.fastq.gz'} in2={mate2+'.fastq.gz'} out={mate_id+'_merged.fastq'} outu={mate_id+'_merge_failed.fastq'} > {wd}stdout.txt\"")
        os.system(f"ls -l > stdout.log")


#stdout_logger = logging.getLogger('STDOUT')
#sl = StreamToLogger(stdout_logger, logging.INFO)
#sys.stdout = sl
#
#stderr_logger = logging.getLogger('STDERR')
#sl = StreamToLogger(stderr_logger, logging.ERROR)
#sys.stderr = sl


#path = f"{wd}stdout.txt"
#sys.stdout = open(path, 'w')