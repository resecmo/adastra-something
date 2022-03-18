import json
import requests
from pathlib import Path
from subprocess import run
from urllib.parse import quote_plus

import main
import bedwriter


def download_adastra_snps():
    r = requests.get(f'https://adastra.autosome.ru/api/v4/browse/tf')
    if r.status_code == 200:
        resp = json.loads(r.text)
        genes = resp['results']
    else:
        raise Exception("Unable to browse tf, got status code {r.status_code}")

    Path("adastra_snps").mkdir(exist_ok=True)
    files = []
    for gene in genes:
        gene_name = gene['gene_name']
        try:
            positions = bedwriter.fetch_adastra_snps(quote_plus(gene_name))
        except Exception as e:  # TODO specify exception
            print(f"got exception in gene {gene_name}, skipping:\n{e}")
            continue
        bed_filename = f"adastra_snps/{gene_name}.bed"
        main.write_positions_in_bed(positions, bed_filename)
        files.append((gene_name, bed_filename))
        print(f"{gene_name} SNP positions ({len(positions)} entries) fetched and written in {bed_filename}")
    return files


def intersect(files):
    result = []
    for gene, bed in files:
        with open("config.json") as config_file:
            config = json.load(config_file)
        cistrome_bed = config["bed_dir"] + f"/{config['gene_name']}_HUMAN.{config['quality']}.bed"
        intersections = run(["bedtools", "intersect",  "-a", cistrome_bed, "-b", bed], capture_output=True)
        intersections = intersections.stdout.decode().split('\n')[:-1]
        intersections = [tuple(intersec.split('\t')) for intersec in intersections]
        result.append((gene, intersections))
    return result


if __name__ == "__main__":
    with open('intersections.txt', 'w') as f:
        f.write(str(intersect(download_adastra_snps())))

