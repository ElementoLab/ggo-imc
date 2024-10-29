import yaml
from pathlib import Path

def load(filename: Path):
    with open(filename, "r") as stream:
        try:
            metadata = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    return metadata