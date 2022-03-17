import json
import requests
from os import system


with open("config.json") as config_file:
    config = json.load(config_file)


def fetch_all_tfs():
    r = requests.get(config['base_api_url'] + f'browse/tf')
    if r.status_code == 200:
        resp = json.loads(r.text)
    else:
        #TODO choose or create proper exception
        raise Exception(f"Tried to browse TFs, got status code: {r.status_code}")
    tfs = resp["results"]
    return [entry['name'] for entry in tfs]


def fetch_adastra_snps(gene_name):
    r = requests.get(config['base_api_url'] + 
                     f'search/snps/advanced?page=1&size=1&offset=0&transcription_factors={gene_name}')
    if r.status_code == 200:
        resp = json.loads(r.text)
        n_snps = resp['total']
    else:
        #TODO choose or create proper exception
        raise Exception(f"Searched for {gene_name}, got status code: {r.status_code}")
    
    #print(n_snps)
    positions = []
    for offset in range(0, n_snps, 100):  #TODO widen
        url = config['base_api_url'] + \
              f'search/snps/advanced?page=1&size=100&offset={offset}&transcription_factors={gene_name}'
        result = requests.get(url)
        snps_json = json.loads(result.text)['results']
        for snp in snps_json:
            positions.append((snp['chromosome'], snp['position']))
    print(positions)
    
    return positions


def write_tuples_in_bed(tuples, filename, append=False):
    file_mode = 'a' if append else 'w'
    with open(filename, file_mode) as bed:
        for tupl in tuples:
            print(*tupl, file=bed, sep="\t")


def write_positions_in_bed(positions, filename, label=None, append=False):
    # chr is pos[0], loc is pos[1]
    if label is None:
        write_tuples_in_bed(map(lambda pos: (pos[0], pos[1]-1, pos[1]), positions),
                            filename, append)
    else:
        write_tuples_in_bed(map(lambda pos: (pos[0], pos[1]-1, pos[1], label), positions),
                            filename, append)


def intersect_adastra_bed(gene_name, path_to_bed):
    adastra_snps = fetch_adastra_snps(gene_name)
    print(f"{len(adastra_snps)} SNPs in ADASTRA")

    adastra_bed_filename = f"adastra-{gene_name}.bed"
    write_positions_in_bed(adastra_snps, adastra_bed_filename)
    print(f"Written SNPs in {adastra_bed_filename}")
    print("Intersections:")

    system(f"bedtools intersect -a {path_to_bed} -b {adastra_bed_filename}")


if __name__ == "__main__":
    #with open("config.json") as config_file:
    #    config = json.load(config_file)
    path_to_bed = config["bed_dir"] + f"/{config['gene_name']}.{config['quality']}.bed"
    intersect_adastra_bed(config["gene_name"], path_to_bed)

