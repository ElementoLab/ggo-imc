from tqdm import tqdm
import requests

def download(url = 'https://wcm.box.com/shared/static/mdntp2xf9tjobxeidkw93mg8jysb7nh9.yml', path = 'metadata/gg_config.yml')

    response = requests.get(url, stream=True)

    with open(path, "wb") as handle:
        for data in tqdm(response.iter_content()):
            handle.write(data)

download()