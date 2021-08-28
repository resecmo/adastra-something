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
adastra_snps = fetch_adastra_snps(config["gene_name"])
print(f"{len(adastra_snps)} SNPs in ADASTRA")
print(adastra_snps)


with open(config["path_to_bed"]) as f:
    bed_lines = f.readlines()
    print(f"{len(bed_lines)} entries in BED")
    for bed_line in bed_lines:
        l = bed_line.split("\t")
        l[-1] = l[-1][:-1]
        for snp in adastra_snps:
            if (snp[0] == l[0]) and (int(l[1]) <= snp[1] <= int(l[2])):
                print(snp, l)


