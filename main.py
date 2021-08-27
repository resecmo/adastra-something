import json
import requests

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
print(fetch_adastra_snps(config["gene_name"]))


with open(config["path_to_bed"]) as f:
    bed_lines = f.readlines()


