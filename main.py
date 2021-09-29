import json
import requests
from os import system


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


def write_positions_in_bed(positions, filename):
    with open(filename, "w") as bed:
        for chrom, loc in positions:
            print(chrom, loc-1, loc, file=bed, sep="\t")


def intersect_adastra_bed(gene_name, path_to_bed):
    adastra_snps = fetch_adastra_snps(gene_name)
    print(f"{len(adastra_snps)} SNPs in ADASTRA")

    adastra_bed_filename = f"adastra-{gene_name}.bed"
    write_positions_in_bed(adastra_snps, adastra_bed_filename)
    print(f"Written SNPs in {adastra_bed_filename}")
    print("Intersections:")

    system(f"bedtools intersect -a {path_to_bed} -b {adastra_bed_filename}")


if __name__ == "__main__":
    with open("config.json") as config_file:
        config = json.load(config_file)
    path_to_bed = config["bed_dir"] + f"/{config['gene_name']}_HUMAN.{config['quality']}.bed"
    intersect_adastra_bed(config["gene_name"], path_to_bed)

