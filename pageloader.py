import json
import requests


def PageLoader(url: str, size: int = 50):
    """A generator that splits ADASTRA API query into chunks and
    yields results of these sub-requests in JSON.
    If any response has status code different from 200, an exception
    is raised.
    """
    if not (1 <= size <= 1000):
        raise ValueError("size should be in range [1;1000]")

    r = requests.get(url + '?size=1')
    if r.status_code != 200:
        raise Exception(f"Got HTTP status code {r.status_code}")
    total_size = json.loads(r.text)['total']

    for offset in range(0, total_size, 50):
        page_size = min(size, total_size-offset)
        r = requests.get(url + f'?offset={offset}&size={page_size}')
        if r.status_code != 200:
            raise Exception(f"Got HTTP status code {r.status_code}")
        yield json.loads(r.text)
