import json 

with open("config.json") as config:
    print(json.load(config))

