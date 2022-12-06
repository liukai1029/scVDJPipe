from os import listdir
from os.path import isfile, join
import os
import logging


#wd = "/export/home/kliu6/projects/single_B_miseq/FASTQ_subset_test/"
#wd = "/export/home/kliu6/projects/single_B_miseq/FASTQ/"
wd = "/export/home/kliu6/projects/single_B_miseq/annotate_test/"





#logging.basicConfig(filename="/fcrbiouatappn01/home/kliu6/projects/single_B_miseq/FASTQ_subset_test/log.txt",
logging.basicConfig(filename=f"{wd}log.txt",
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)

def call_pyir(f):
        f_root = f.replace("_merged.fastq", "")
        zipped_name = f_root+".tsv.gz"
        unzipped_name = f_root+".tsv"
        os.system(f"/bin/bash -c \"pyir -t fastq --igdata /pyir_data  -r Ig -s mouse -x /pyir/pyir/data/bin/igblastn_linux -o {f_root} --outfmt tsv {f} \"")
        os.system(f"gzip -cd {zipped_name} > {unzipped_name}")

        
fl = [f for f in listdir(wd) if isfile(join(wd, f))]
fl = [f for f in fl if f[-13:] == "_merged.fastq" ]
fl = sorted(fl)
fp = [wd+f for f in fl]


logging.info(f"all merged files identified: {fp}")


while fp:
    input = fp.pop(0)
    logging.info(f"start pyir for: {input}")
    call_pyir(input)    
    logging.info(f"complete pyir for: {input}")
