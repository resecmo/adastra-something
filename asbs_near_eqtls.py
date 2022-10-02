import glob
import subprocess
from time import strftime, gmtime
import typing

import pandas as pd
import requests
# from urllib3.util import Retry

import bedwriter
from pageloader import PageLoader


def download_eqtls(gene_name: str):
    eqtl_tuples = []
    adastra_url = f'https://adastra.autosome.org/api/v5/search/snps/eqtl_gene_name/{gene_name}'
    for eqtls_json in PageLoader(adastra_url, size=500):
        eqtl_tuples.extend((eqtl['chromosome'], eqtl['position']-1, eqtl['position'])
                           for eqtl in eqtls_json['results'])
    bedwriter.write_tuples_in_bed(eqtl_tuples, f"eqtls/{gene_name}-eqtl.bed")  # write_positions?


def download_eqtls_multigene(gene_names: typing.Iterable[str]):
    s = requests.Session()
    # retries = Retry(total=3, status_forcelist=[500])
    # a = requests.adapters.HTTPAdapter(max_retries=retries)
    # s.mount('https://', a)

    for i, gene_name in enumerate(gene_names):
        eqtl_tuples = []
        adastra_url = f'https://adastra.autosome.org/api/v5/search/snps/eqtl_gene_name/{gene_name}'
        for eqtls_json in PageLoader(adastra_url, size=500, session=s):
            eqtl_tuples.extend((eqtl['chromosome'], eqtl['position']-1, eqtl['position'])
                               for eqtl in eqtls_json['results'])
        bedwriter.write_tuples_in_bed(eqtl_tuples, f"eqtls/{gene_name}-eqtl.bed")  # write_positions?

        if (i+1) % 10 == 0:
            print(f'{strftime("%X %x %Z", gmtime())}: '
                  f'eQTL data has been downloaded for {i+1} out of {len(all_genes)} genes')


def find_asbs_near_eqtls(gene_name, window_halfwidth: int = 2000, asbs_scope="all"):  # todo rename
    """
    Searches for ASBs in vicinity of eQTLs with target `gene_name`; 
    possibly these TFs regulate `gene_name` expression?
    
    `asb_scope` defines range of search for asbs:
    `"all"` for all TFs; 
    `"crohn"`, `"ulcerative_colitis"`, `"type_ii_diabetes"`
    or `"breast_cancer"` for TFs from the article.
    """

    if (asbs_scope is None) or (asbs_scope == "all"):
        asbs_filename = "adastra_snps/all_tfs.bed"
    elif asbs_scope in ["crohn", "ulcerative_colitis", "type_ii_diabetes", "breast_cancer"]:
        asbs_filename = f"adastra_snps/{asbs_scope}_tfs.bed"
    else:
        raise ValueError("asbs_scope should be one of: \"all\", \"diabetes\"")

    # widen sort and merge intersecting intervals
    widened_bed = subprocess.run(["bedtools", "slop", "-b", f'{window_halfwidth}',
                                  "-i", f"eqtls/{gene_name}-eqtl.bed",
                                  "-g", "/usr/share/bedtools/genomes/human.hg38.genome"],
                                 capture_output=True).stdout
    with open(f"eqtls/window_beds/{gene_name}-eqtl.bed", 'wb') as f:
        f.write(widened_bed)

    sorted_bed = subprocess.run(["bedtools", "sort", "-i", f"eqtls/window_beds/{gene_name}-eqtl.bed"],
                                capture_output=True).stdout
    with open(f"eqtls/window_beds/{gene_name}-eqtl.bed", 'wb') as f:
        f.write(sorted_bed)

    merged_bed = subprocess.run(["bedtools", "merge", "-i", f"eqtls/window_beds/{gene_name}-eqtl.bed"],
                                capture_output=True).stdout
    with open(f"eqtls/window_beds/{gene_name}-eqtl-merged.bed", 'wb') as f:
        f.write(merged_bed)

    intersection = subprocess.run(["bedtools", "intersect", "-a", asbs_filename,
                                   "-b", f"eqtls/window_beds/{gene_name}-eqtl-merged.bed", "-wa"],
                                  capture_output=True).stdout

    return f"ASBs near eQTLs with target {gene_name}:\n{intersection.decode()}"


def prepare_bed_for_tfs_subset(trait_name: str, gene_names_list, p_thr):
    """
    Selects all statistically significant ASBs of provided genes
    from a local copy of ADASTRA and saves them to a .bed file
    f"adastra_snps/{trait_name}_tfs.bed"
    """

    # Filename of the created .bed
    bed_filename = f"adastra_snps/{trait_name}_tfs.bed"
    with open(bed_filename, 'w'):
        pass
    for gene_name in gene_names_list:
        for gene_bed_filename in glob.glob(f"/home/resecmo/vigg/release_BillCipher/release_dump/TF/*{gene_name}_*"):
            df = pd.read_csv(gene_bed_filename, sep='\t', usecols=[0, 1, 2, 14, 15])  # 2 is id
            snps = []
            for i in range(df.shape[0]):
                if (df["fdrp_bh_ref"][i] < p_thr) or (df["fdrp_bh_alt"][i] < p_thr):
                    snps.append((df["#chr"][i], df["pos"][i] - 1, df["pos"][i], df["ID"][i], f"{gene_name}_HUMAN"))
            bedwriter.write_tuples_in_bed(snps, bed_filename, append=True)


if __name__ == '__main__':

    # genes list is from https://doi.org/10.1101/2021.06.08.21258515. p.17
    genes = {
        'ldl': ['APOB', 'APOC2', 'APOE', 'LDLR', 'LPL', 'PCSK9'],
        'hdl': ['ABCA1', 'APOA1', 'CETB', 'LIPC', 'LIPG', 'PLTP', 'SCARB1'],
        'height': [
           'ANTXR1', 'ATR', 'BLM', 'CDC6', 'CDT1', 'CENPJ', 'COL1A1', 'COL1A2', 'COMP', 'CREBBP', 'DNA2',
           'DTDST', 'EP300', 'EVC', 'EVC2', 'FBN1', 'FGFR3', 'FKBP10', 'GHR', 'KRAS', 'NBN', 'NIPBL', 'ORC1',
           'ORC4', 'ORC6L', 'PCNT', 'PLOD2', 'PTPN11', 'RAD21', 'RAF1', 'RECQL4', 'RIT1',
           'RNU4ATAC should remove snRNA', 'ROR2', 'SLC26A2', 'SMAD4', 'SMC3 milder form of trait, remove',
           'SOS1 same as above', 'SRCAP', 'WRN'],
        'pressure': ['KCNJ1', 'SLC12A1', 'SLC12A3', 'WNK1', 'WNK4'],
        'crohn': [  # K50
            'ATG16L1', 'CARD9', 'IL10', 'IL10RA', 'IL10RB', 'IL23R', 'IRGM', 'NOD2', 'PRDM1', 'PTPN22', 'RNF'],
        'ulcerative_colitis': ['ATG16L1', 'CARD9', 'IL23R', 'IRGM', 'PRDM1', 'PTPN22', 'RNF186'],  # K51
        'type_ii_diabetes': [  # E11
            'ABCC8', 'BLK', 'CEL', 'EIF2AK3', 'GATA4', 'GATA6', 'GCK', 'GLIS3', 'HNF1A', 'HNF1B',
            'HNF4A', 'IER3IP1', 'INS', 'KCNJ11', 'KLF11', 'LMNA', 'NEUROD1', 'NEUROG3', 'PAX4', 'PDX1',
            'PPARG', 'PTF1A', 'RFX6', 'SLC19A2', 'SLC2A2', 'WFS1', 'ZFP57'],
        'breast_cancer': [  # C50? Malignant neoplasms of breast
            'AKT1', 'ARID1A', 'ATM', 'BRCA1', 'BRCA2', 'CBFB', 'CDH1',
            'CDKN1B', 'CHEK2', 'CTCF', 'ERBB2', 'ESR1', 'FGFR2', 'FOXA1',
            'GATA3', 'GPS2', 'HS6ST1', 'KMT2C', 'KRAS', 'LRRC37A3', 'MAP2K4',
            'MAP3K1', 'NCOR1', 'NF1', 'NUP93', 'PALB2', 'PIK3CA', 'PTEN',
            'RB1', 'RUNX1', 'SF3B1', 'STK11', 'TBX3', 'TP53', 'ZFP36L1']
    }
    full_trait_names = {'ldl': 'LDL', 'hdl': 'HDL', 'height': 'Height',
                        'pressure': 'Blood pressure (systolic and diastolic)',
                        'crohn': 'Crohn disease',
                        'ulcerative_colitis': 'Ulcerative colitis',
                        'type_ii_diabetes': 'Type II diabetes',
                        'breast_cancer': 'Breast cancer (selected using MutPanning26)'}
    chosen_traits = ['ldl', 'hdl', 'pressure', 'crohn', 'ulcerative_colitis',  # without height - TODO broken genes
                     'type_ii_diabetes', 'breast_cancer']
    all_genes = set.union(*(set(genes_sublist) for trait, genes_sublist in genes.items() if trait in chosen_traits))

    download_eqtls_multigene(all_genes)
    print('eQTL data has been downloaded')

    # P-value threshold for ADASTRA's ASBs
    p_thr = 0.05
    for trait in ['crohn', 'ulcerative_colitis', 'type_ii_diabetes', 'breast_cancer']:
        trait_genes = genes[trait]
        prepare_bed_for_tfs_subset(trait, trait_genes, p_thr)

    window_sizes = [500 * 2**i for i in range(16)]
    for trait in ['crohn', 'ulcerative_colitis', 'type_ii_diabetes', 'breast_cancer']:
        for wsize in window_sizes:
            result = ''
            for gene_name in genes[trait]:
                result += find_asbs_near_eqtls(gene_name, window_halfwidth=wsize, asbs_scope=trait)
            with open(f'asbs_near_eqtls/{trait}/w{wsize}_p{p_thr}.txt', 'w') as f:
                f.write(result)
            print(f'{wsize} done at {strftime("%X %x %Z", gmtime())}')
        
