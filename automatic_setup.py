
import os
from git import Repo
import shutil
import gzip
import zipfile
import requests
import math
from pathlib import Path
from multiprocessing import Pool, cpu_count
from tqdm import tqdm


def unzip(path):
    goal = "/".join(path.split(".zip")[0].split("/")[:-1])
    with zipfile.ZipFile(path, "r") as f:
        f.extractall(goal)


def ungunzip(path):
    goal = path.split(".gz")[0]
    with gzip.open(path, 'rb') as f_in, open(goal, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)


def downloader(task, multithread=True):
    path, url = task
    file_path = "{path}/{filename}".format(
        path=path, filename=url.split('/')[-1])
    if not os.path.exists(file_path) and not os.path.exists(".".join(
            file_path.split(".")[:-1])):
        r = requests.get(url, stream=True)
        total_size = int(r.headers.get('content-length', 0))
        block_size = 1024**2
        wrote = 0
        with open(file_path, 'wb') as f:
            for data in tqdm(
                    r.iter_content(block_size),
                    disable=multithread,
                    total=math.ceil(total_size // block_size),
                    unit='MB',
                    leave=False,
                    unit_scale=True):
                wrote = wrote + len(data)
                f.write(data)
        if file_path.endswith(".zip"):
            unzip(file_path)
            os.remove(file_path)
        elif file_path.endswith(".gz"):
            ungunzip(file_path)
            os.remove(file_path)


url = "https://github.com/gagneurlab/Manuscript_Avsec_Bioinformatics_2017"
path = url.split("/")[-1]


pip_requirements = "{path}/python3_requirements.txt".format(path=path)
r_installer = "{path}/R_installer.r".format(path=path)
raw_path = "{path}/data/eclip/raw".format(path=path)
fasta_path = "{path}/fasta".format(path=path)
simdna_url = "https://github.com/kundajelab/simdna/archive/0.2.zip"
simdna_path = "./simdna-0.2"
metadata_path = "{raw_path}/metadata.tsv".format(raw_path=raw_path)
metadata_url = "https://github.com/gagneurlab/Manuscript_Avsec_Bioinformatics_2017/files/2447032/metadata.tsv.zip"
fasta_url = "https://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
annotations_url = "https://github.com/gagneurlab/Manuscript_Avsec_Bioinformatics_2017/files/2484555/gencode.v25.annotation.gtf.zip"
target_fasta = "{fasta_path}/GRCh38.p7.genome.fa".format(fasta_path=fasta_path)
file_source = "{path}/Scripts/RBP/Eclip/files.txt".format(path=path)


print("Cloning repository.")
if not os.path.exists(path): 
    Repo.clone_from(url, path)


print("Requiring simdna==0.2, as it is no longer available on pipy.")
if not os.path.exists(simdna_path): 
    downloader(("./", simdna_url))
print("Please run from a terminal in this directory the following:\n")
print("\033[1m  pip install -e {simdna_path}  \033[0m".format(simdna_path=simdna_path))
input("When done, press any key to continue... ")
shutil.rmtree(simdna_path)


print("Now we proceed to install the required python packages.")
print("Please run from a terminal in this directory the following:\n")
print("\033[1m  pip install -r {pip_requirements}  \033[0m".format(pip_requirements=pip_requirements))
input("When done, press any key to continue... ")


print("Now we proceed to install the required R packages.")
print("Please run from a terminal in this directory the following:\n")
print("\033[1m Rscript {r_installer} \033[0m".format(r_installer=r_installer))
input("When done, press any key to continue... ")


print("Downloading tab files.")
with open(file_source, "r") as f:
    urls = f.read().split("\n")
    
Path(raw_path).mkdir(parents=True, exist_ok=True)

jobs = [(raw_path, url) for url in urls]

with Pool(cpu_count()) as p:
    list(tqdm(p.imap(downloader, jobs), total=len(jobs)))


print("Switching metadata.tsv file.")
if os.path.exists(metadata_path):
    os.remove(metadata_path)
downloader((raw_path, metadata_url), multithread=False)


print("Downloading fasta.")
Path(fasta_path).mkdir(parents=True, exist_ok=True)
downloader((fasta_path, fasta_url), multithread=False)


shutil.move("{fasta_path}/{expected_fasta_name}".format(
    fasta_path=fasta_path,
    expected_fasta_name=fasta_url.split("/")[-1]
), target_fasta)


print("Downloading annotations.")
downloader((fasta_path, annotations_url), multithread=False)


print("You should now be able to run snakemake.")

