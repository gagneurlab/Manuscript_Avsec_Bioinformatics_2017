import os
import pandas as pd
import numpy as np
from tqdm import tqdm_notebook as tqdm


def clear_file_chromosomes(filename: str, target: str):
    """Removes rows with chromosomes beginning with given target.
        filename:str, the path to the file to clear.
        target:str, the target for the rows to remove.
    """
    df = pd.read_csv(filename, sep='\t', header=None)
    a = np.array(df)
    mask = np.array([row.startswith(target) for row in a[:, 0]])
    pd.DataFrame(a[~mask]).to_csv(filename, sep='\t', header=None)
    

def clear_directory_chromosomes(directory: str, target: str = "chrUn_"):
    """Removes rows with chromosomes beginning with given target.
        directory:str, the path to the directory to clear.
        target:str, the target for the rows to remove.
    """
    [
        clear_file_chromosomes(
            "{directory}/{file}".format(directory=directory, file=file),
            target) for file in tqdm(next(os.walk(directory))[2]) if file.endswith(".tab")
    ]
    
if __name__ == '__main__':
    clear_directory_chromosomes("../../../data/eclip/raw")
