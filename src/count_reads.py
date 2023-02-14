from os import listdir
from os.path import isfile, join
import os
import logging
import pandas as pd
import numpy as np

#wd = "/export/home/kliu6/projects/single_B_miseq/FASTQ_subset_test/"
#wd = "/export/home/kliu6/projects/single_B_miseq/FASTQ/"
#wd = "/export/home/kliu6/projects/single_B_miseq/rc_test/"
#wd = "/export/home/kliu6/projects/single_B_miseq/FASTQ_human/"
wd = "/export/home/kliu6/projects/single_B_miseq/FC1A/"





#logging.basicConfig(filename="/fcrbiouatappn01/home/kliu6/projects/single_B_miseq/FASTQ_subset_test/log.txt",
logging.basicConfig(filename=f"{wd}log_read_count.txt",
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)

required_col = ['sequence_id',
 'sequence',
 'locus',
 'stop_codon',
 'vj_in_frame',
 'productive',
 'v_call',
 'd_call',
 'j_call',
 'sequence_alignment',
 'germline_alignment',
 'sequence_alignment_aa',
 'germline_alignment_aa',
 'fwr1',
 'fwr1_aa',
 'cdr1',
 'cdr1_aa',
 'fwr2',
 'fwr2_aa',
 'cdr2',
 'cdr2_aa',
 'fwr3',
 'fwr3_aa',
 'fwr4',
 'fwr4_aa',
 'cdr3',
 'cdr3_aa',
 'v_family',
 'd_family',
 'j_family',
 'cdr3_aa_length']

ordered_col = ['count',
 'frequency',
 'sequence',
 'locus',
 'stop_codon',
 'vj_in_frame',
 'productive',
 'v_call',
 'd_call',
 'j_call',
 'sequence_alignment',
 'germline_alignment',
 'sequence_alignment_aa',
 'germline_alignment_aa',
 'fwr1',
 'fwr1_aa',
 'cdr1',
 'cdr1_aa',
 'fwr2',
 'fwr2_aa',
 'cdr2',
 'cdr2_aa',
 'fwr3',
 'fwr3_aa',
 'fwr4',
 'fwr4_aa',
 'cdr3',
 'cdr3_aa',
 'v_family',
 'd_family',
 'j_family',
 'cdr3_aa_length'
]

def count_reads(tsv_fp):
    df = pd.read_csv(tsv_fp, sep="\t")
    keep = df.copy()
    keep = keep[keep["stop_codon"]=="F"]
    keep = keep[keep["productive"]=="T"]
    keep = keep[~keep["v_call"].isna()]
    keep = keep[~keep["j_call"].isna()]
    keep = keep[~keep['fwr1_aa'].isna()]
    keep = keep[~keep['fwr2_aa'].isna()]
    keep = keep[~keep['fwr3_aa'].isna()]
    keep = keep[~keep['fwr4_aa'].isna()]
    keep = keep[~keep['cdr1_aa'].isna()]
    keep = keep[~keep['cdr2_aa'].isna()]
    keep = keep[~keep['cdr3_aa'].isna()]
    action = ['nunique'] + [lambda x:max(set(list(x)), key = list(x).count)]*(len(required_col)-1)
    vdj_keep = keep.groupby('sequence_alignment_aa').agg(dict(zip(required_col, action)))
    new_required_col = ['count' if ele=='sequence_id' else ele for ele in required_col]
    vdj_keep.columns = new_required_col
    vdj_keep.sort_values(by='count', ascending=False, inplace=True)
    vdj_keep = vdj_keep[['*' not in ele for ele in vdj_keep['sequence_alignment_aa'].to_list()]]
    total_read_count = vdj_keep['count'].sum()
    vdj_keep['frequency']=vdj_keep['count']/total_read_count
    vdj_keep = vdj_keep[ordered_col]
    vdj_keep.shape
    #vdj_keep.reset_index(inplace=True)
    fname = tsv_fp.split("/")[-1:][0]
    sname = fname.replace('.tsv','')
    vdj_keep.index = [sname+'_'+str(ele) for ele in list(np.arange(1, len(vdj_keep) + 1))]
    rc_fp = tsv_fp.replace(".tsv","_rc.csv")
    vdj_keep.to_csv(rc_fp)



        
fl = [f for f in listdir(wd) if isfile(join(wd, f))]
fl = [f for f in fl if f[-4:] == ".tsv" ]
fl = sorted(fl)
fp = [wd+f for f in fl]
logging.info(f"all merged files identified: {fp}")

rc_fl = [f for f in listdir(wd) if isfile(join(wd, f))]
rc_fl = [f for f in rc_fl if f[-7:] == "_rc.csv" ]
rc_fl  = sorted(rc_fl )
rc_fp = [wd+f for f in rc_fl]

while fp:
    input = fp.pop(0)
    if input.replace(".tsv","_rc.csv") in rc_fp:
        logging.info(f"Read count file already exist for: {input}")
    else:
        try:
            logging.info(f"start counting reads for: {input}")
            count_reads(input)    
            logging.info(f"complete counting reads for: {input}")
        except KeyError:
            logging.info(f"\n-\n-\n-\n-failed counting reads for: {input}\n-\n-\n-\n-")
            continue






#if "IGH" in df["locus"].to_list():
#    df_h = df[df["locus"]=="IGH"]
#if "IGK" in df["locus"].to_list():
#    df_k = df[df["locus"]=="IGK"]
#if "IGL" in df["locus"].to_list():
#    df_l = df[df["locus"]=="IGL"]







#while fp:
#    input = fp.pop(0)
#    logging.info(f"start pyir for: {input}")
#    call_pyir(input)    
#    logging.info(f"complete pyir for: {input}")
