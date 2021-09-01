import json
import requests
from os import system

with open("config.json") as config_file:
    config = json.load(config_file)


def fetch_adastra_snps(gene_name):
    r = requests.get(f'https://adastra.autosome.ru/api/v3/search/snps/gene_name/{gene_name}')
    if r.status_code == 200:
        resp = json.loads(r.text)
    else:
        print(f"status code: {r.status_code}")
    adastra_snps = resp["results"]
    positions = []
    for snp in adastra_snps:
        #print(json.dumps(snp, indent=2))
        #print(f'{snp["chromosome"]}:{snp["position"]}')
        positions.append((snp["chromosome"], snp["position"]))
    return positions
adastra_snps = fetch_adastra_snps(config["gene_name"])
print(f"{len(adastra_snps)} SNPs in ADASTRA")
print(adastra_snps)
with open(f"adastra-{config['gene_name']}.bed", "w") as out:
    for chrom, loc in adastra_snps:
        print(chrom, loc, loc, file=out, sep="\t")
system("bedtools intersect -a adastra-HNF4A.bed -b ../cistrome_hg38/hg38_cistrome/HNF4A_HUMAN.A.bed")


