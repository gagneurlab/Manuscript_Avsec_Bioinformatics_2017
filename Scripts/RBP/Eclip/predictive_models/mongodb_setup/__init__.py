import json
import os
path = os.path.dirname(os.path.realpath(__file__))

with open("{path}/setup.json".format(path=path), "r") as f:
    data = json.load(f)

host = data["host"]
port = data["port"]

__all__ = ["host", "port"]
