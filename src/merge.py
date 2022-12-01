
from os import listdir
from os.path import isfile, join
import os


wd = "/fcrbiouatappn01/home/kliu6/projects/single_B_miseq/FASTQ/"

fl = [f for f in listdir(wd) if isfile(join(wd, f))]

fl = sorted(fl)

test_f = fl[40:42]

test_fp = [wd+f for f in test_f]

in_fname = test_fp[0]

ext_name = ".fastq.gz"
mate1 = in_fname.replace(ext_name, "")

mate2 = mate1.replace("_R1_","_R2_")



if mate2+ext_name in test_fp:
    print("yes")
    mate_id = mate1.replace("_R1_001", "")
    print(f"/bin/bash -c \"/fcrbiouatappn01/home/kliu6/packages/bbmap/bbmerge.sh in1={mate1+'.fastq.gz'} in2={mate2+'.fastq.gz'} out={mate_id+'_merged.fastq'} outu={mate_id+'_merge_failed.fastq'}\"")
    os.system(f"/bin/bash -c \"/fcrbiouatappn01/home/kliu6/packages/bbmap/bbmerge.sh in1={mate1+'.fastq.gz'} in2={mate2+'.fastq.gz'} out={mate_id+'_merged.fastq'} outu={mate_id+'_merge_failed.fastq'}\"")


