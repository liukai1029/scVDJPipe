from os import listdir
from os.path import isfile, join, exists
import os
import logging
import pandas as pd
import numpy as np

#wd = "/export/home/kliu6/projects/single_B_miseq/FASTQ_subset_test/"
wd = "/export/home/kliu6/projects/single_B_miseq/FASTQ/"
#wd = "/export/home/kliu6/projects/single_B_miseq/rc_test/"





#logging.basicConfig(filename="/fcrbiouatappn01/home/kliu6/projects/single_B_miseq/FASTQ_subset_test/log.txt",
logging.basicConfig(filename=f"{wd}log_agg.txt",
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)

temp_fp = wd+"main_heavy_light_chains_temp.csv"
agg_fp = wd+"Main_Heavy_Light_Chains_Table.csv"




def choose_major_chain(rc_fp, agg_fp):
    rc = pd.read_csv(rc_fp, index_col = 0)
    rc_major = rc.copy()
    rc_major = rc_major[rc_major['frequency']>0.01]

    rc_h = pd.DataFrame(columns=required_col)
    rc_k = pd.DataFrame(columns=required_col)
    rc_l = pd.DataFrame(columns=required_col)

    if "IGH" in rc_major["locus"].to_list():
        rc_h = rc_major[rc_major["locus"]=="IGH"]
        h_criteria = rc_h["frequency"] > rc_h.iloc[0]["frequency"]/5
        rc_h = rc_h[h_criteria]

    if "IGK" in rc_major["locus"].to_list():
        rc_k = rc_major[rc_major["locus"]=="IGK"]
        k_criteria = rc_k["frequency"] > rc_k.iloc[0]["frequency"]/5
        rc_k = rc_k[k_criteria]

    if "IGL" in rc_major["locus"].to_list():
        rc_l = rc_major[rc_major["locus"]=="IGL"]
        l_criteria = rc_l["frequency"] > rc_l.iloc[0]["frequency"]/5
        rc_l = rc_l[l_criteria]

    try:
        major_chains = pd.concat([rc_h, rc_k, rc_l])
    except TypeError:
        import pdb; pdb.set_trace()

    # replace "," with "/" in vdj calls with multiple match
    major_chains['v_call'] = [str(ele).replace(",","/") for ele in  major_chains['v_call'].to_list()]
    major_chains['d_call'] = [str(ele).replace(",","/") for ele in  major_chains['d_call'].to_list()]
    major_chains['j_call'] = [str(ele).replace(",","/") for ele in  major_chains['j_call'].to_list()]

    return major_chains






# identify all tsv pyir files
fl = [f for f in listdir(wd) if isfile(join(wd, f))]
fl = [f for f in fl if f[-4:] == ".tsv" ]
fl = sorted(fl)
fp = [wd+f for f in fl]
logging.info(f"all merged files identified: {fp}")

# identify all read count files
rc_fl = [f for f in listdir(wd) if isfile(join(wd, f))]
rc_fl = [f for f in rc_fl if f[-7:] == "_rc.csv" ]
rc_fl  = sorted(rc_fl )
rc_fp = [wd+f for f in rc_fl]
logging.info(f"all merged files identified: {rc_fp}")

sample = pd.read_csv(rc_fp[0], index_col = 0)
required_col = sample.columns.to_list()


if not exists(agg_fp):
    with open(agg_fp, "w") as f:
        f.write('ID,')
        for ele in required_col:
            f.write(str(ele)+',')
        f.write('\n')
else:
    logging.info(f"aggregated major chain file already exist: {agg_fp}")



# loop through all read count files
while rc_fp:
    input = rc_fp.pop(0)
    if input.replace("_rc.csv",".tsv") in fp:
        try:
            logging.info(f"Read count file exist: {input}")
            logging.info(f"Start major chain selection: {input}")
            major_chains_by_input = choose_major_chain(input, agg_fp)

            if exists(agg_fp):
                with open (agg_fp, "a") as f:
                    for i in range(major_chains_by_input.shape[0]):
                        record = major_chains_by_input.iloc[i].to_list()
                        f.write(str(major_chains_by_input.index[i])+',')
                        for ele in record:
                            f.write(str(ele)+',')
                        f.write('\n')

        except KeyError:
            logging.info(f"\n-\n-\n-q\n-failed getting major chain for: {input}\n-\n-\n-\n-")
            continue    
    else:
        logging.info(f"Major chain file not available for: {input.replace('_rc.csv','.tsv')}")






