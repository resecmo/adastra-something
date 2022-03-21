import glob
import json
import requests
import subprocess

import pandas as pd

import bedwriter


def find_asbs_near_eqtls(gene_id, gene_name, asbs_scope="all"):  # todo rename
    """
    Searches for ASBs in vicinity of eQTLs with target `gene_name`; 
    possibly these TFs regulate `gene_name` expression?
    
    `asb_scope` defines range of search for asbs:
    `"all"` for all TFs; 
    `"diabetes"` for TFs from the article.
    """
    
    
    if asbs_scope == "all":
        asbs_filename = "adastra_snps/all_tfs.bed"
    elif asbs_scope == "diabetes":
        asbs_filename = "adastra_snps/diabetes_tfs.bed"
    else:
        raise ValueError("asbs_scope should be one of: \"all\", \"diabetes\"") 
    
    
    request = requests.get(f'https://adastra.autosome.ru/api/v4/search/snps/eqtl_gene_id/' + gene_id)
    eqtls = json.loads(request.text)['results']



    window_halfwidth = 2000
    eqtl_tuples = []
    for eqtl in eqtls:
        eqtl_tuples.append((eqtl['chromosome'],
                            eqtl['position'] - 1 - window_halfwidth,
                            eqtl['position'] + window_halfwidth))
    bedwriter.write_tuples_in_bed(eqtl_tuples, f"eqtls/{gene_name}-eqtl.bed")
    #print(*eqtl_tuples, sep='\n', end='\n\n')

    # sort and merge intersecting intervals
    sorted_bed = subprocess.run(["bedtools", "sort", "-i", f"eqtls/{gene_name}-eqtl.bed"],
                                capture_output=True).stdout
    with open(f"eqtls/{gene_name}-eqtl.bed", 'wb') as f:
        f.write(sorted_bed)

    merged_bed = subprocess.run(["bedtools", "merge", "-i", f"eqtls/{gene_name}-eqtl.bed"],
                                capture_output=True).stdout
    with open(f"eqtls/{gene_name}-eqtl-merged.bed", 'wb') as f:
        f.write(merged_bed)

    #print(merged_bed.decode())


        
    intersection = subprocess.run(["bedtools", "intersect", "-a", asbs_filename,
                                   "-b", f"eqtls/{gene_name}-eqtl-merged.bed", "-wa"], 
                                  capture_output=True).stdout
    
    print(f"ASBs near eQTLs with target {gene_name}:")
    print(intersection.decode())

if __name__ == '__main__':
    #find_asbs_near_eqtls('ENSG00000101076', 'HNF4A') # HNF4A
    #find_asbs_near_eqtls('ENSG00000204644', 'ZFP57')

    # genes list is from https://doi.org/10.1101/2021.06.08.21258515. p.17
    # 
    # todo is GATA4_alt needed?
    genes = \
    '''ABCC8 ENSG00000006071
    BLK ENSG00000136573
    CEL ENSG00000170835
    EIF2AK3 ENSG00000172071
    GATA4 ENSG00000136574 
    GATA4_alt ENSG00000285109
    GATA6 ENSG00000141448
    GCK ENSG00000106633
    GLIS3 ENSG00000107249
    HNF1A ENSG00000135100
    HNF1B ENSG00000275410
    HNF4A ENSG00000101076
    IER3IP1 ENSG00000134049
    INS ENSG00000254647
    KCNJ11 ENSG00000187486
    KLF11 ENSG00000172059
    LMNA ENSG00000160789
    NEUROD1 ENSG00000162992
    NEUROG3 ENSG00000122859
    PAX4 ENSG00000106331
    PDX1 ENSG00000139515
    PPARG ENSG00000132170
    PTF1A ENSG00000168267
    RFX6 ENSG00000185002
    SLC19A2 ENSG00000117479
    SLC2A2 ENSG00000163581
    WFS1 ENSG00000109501
    ZFP57 ENSG00000204644'''
    
    genes = [s.strip().split() for s in genes.split('\n')]
    gene_names = [gene[0] for gene in genes]
    #print(gene_names)
    
    
    # P-value threshold for ADASTRA's ASBs
    p_thr = 0.05
            
    with open("adastra_snps/diabetes_tfs.bed", 'w'):
        pass
    for gene_name in gene_names:
        #print(gene_name)
        #print(glob.glob(f"/home/resecmo/vigg/release_Zanthar/release_dump/TF/*{gene_name}_*"))
        for gene_bed_filename in glob.glob(f"/home/resecmo/vigg/release_Zanthar/release_dump/TF/*{gene_name}_*"):
            df = pd.read_csv(gene_bed_filename, sep='\t', usecols=[0, 1, 2, 14, 15])  # 2 is id
            snps = []
            for i in range(df.shape[0]):
                if ((df["fdrp_bh_ref"][i] < p_thr) or (df["fdrp_bh_alt"][i] < p_thr)):
                    snps.append((df["#chr"][i], df["pos"][i] - 1, df["pos"][i], df["ID"][i], f"{gene_name}_HUMAN"))
            bedwriter.write_tuples_in_bed(snps, "adastra_snps/diabetes_tfs.bed", append=True)



    
    #print(genes)
    for gene_name, gene_id in genes: 
        find_asbs_near_eqtls(gene_id, gene_name, asbs_scope="diabetes")
        
    # python asbs_near_eqtls.py >> diabetes_asbs.txt

    
    