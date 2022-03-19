import json
import requests
import subprocess
import bedwriter


def find_asbs_near_eqtls(gene_id, gene_name):  # todo rename
    #searches for ASBs in vicinity of eQTLs with target `gene_name`; 
    #possibly these TFs regulate `gene_name` expression?
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



    intersection = subprocess.run(["bedtools", "intersect", "-a", "adastra_snps/all_tfs.bed",
                                   "-b", f"eqtls/{gene_name}-eqtl-merged.bed", "-wa"], 
                                  capture_output=True).stdout
    print(intersection.decode())

if __name__ == '__main__':
    #find_asbs_near_eqtls('ENSG00000101076', 'HNF4A') # HNF4A
    find_asbs_near_eqtls('ENSG00000204644', 'ZFP57')


'''
ABCC8 ENSG00000006071
BLK ENSG00000136573
CEL ENSG00000170835
EIF2AK3 ENSG00000172071
GATA4 ENSG00000136574 alt:ENSG00000285109
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
ZFP57 ENSG00000204644
'''