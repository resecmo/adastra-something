import json
import requests
from pathlib import Path
from subprocess import run
from urllib.parse import quote_plus

import main

r = requests.get(f'https://adastra.autosome.ru/api/v4/browse/tf')
if r.status_code == 200:
    resp = json.loads(r.text)
    genes = resp['results']
else:
    raise Exception("Unable to browse tf, got status code {r.status_code}")

Path("adastra_snps").mkdir(exist_ok=True)
files = []
for gene in genes[620:625]:
    gene_name = gene['gene_name']
    try:
        positions = main.fetch_adastra_snps(quote_plus(gene_name))
    except Exception as e:  # TODO specify exception
        print(f"got exception in gene {gene_name}, skipping:\n{e}")
        continue
    bed_filename = f"adastra_snps/{gene_name}.bed"
    main.write_positions_in_bed(positions, bed_filename)
    files.append((gene_name, bed_filename))
    print(f"{gene_name} SNP positions ({len(positions)} entries) fetched and written in {bed_filename}")


result = []
for gene, bed in files:
    with open("config.json") as config_file:
        config = json.load(config_file)
    cistrome_bed = config["bed_dir"] + f"/{config['gene_name']}_HUMAN.{config['quality']}.bed"
    print(gene)
    #system(f"bedtools intersect -a {cistrome_bed} -b {bed}")
    intersections = run(["bedtools", "intersect",  "-a", cistrome_bed, "-b", bed], capture_output=True)
    result.append((gene, intersections.stdout))
print(result)

