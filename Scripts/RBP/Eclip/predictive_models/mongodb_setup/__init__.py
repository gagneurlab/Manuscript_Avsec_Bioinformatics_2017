import json

with open("setup.json", "r") as f:
    data = json.load(f)

host = data["host"]
port = data["port"]

__all__ = ["host", "port"]
