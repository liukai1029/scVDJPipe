import os

f_name_list = ["~/projects/single_B_miseq/FASTQ/SAM003330_LIB003809_2_S2_L001_merged.fastq"]


def call_pyir_group(f_name_list):
    for in_fname in f_name_list:
        ext_name = ".fastq"
        core_name = in_fname.replace(ext_name, "")
        zipped_name = core_name+".tsv.gz"
        unzipped_name = core_name+".tsv"
        os.system(f"/bin/bash -c \"pyir -t fastq --igdata /pyir_data  -r Ig -s mouse -x /pyir/pyir/data/bin/igblastn_linux -o {core_name} --outfmt tsv {in_fname} \"")
        #os.system(f"gzip -cd {zipped_name} > {unzipped_name}")


call_pyir_group(f_name_list)    
