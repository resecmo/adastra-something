import json
import requests

with open("config.json") as config_file:
    config = json.load(config_file)
gene_name = config["gene_name"]


r = requests.get(f'https://adastra.autosome.ru/api/v3/search/snps/gene_name/{gene_name}')
if r.status_code == 200:
    resp = json.loads(r.text)
else:
    print(f"status code: {r.status_code}")



print(json.dumps(resp, indent=2))

