import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import sys
import subprocess
import os
from tqdm import tqdm


if __name__ == '__main__':
    
    args = sys.argv
    datafile = args[1]
    sigmafile = args[2]
    xifile = args[3]
    threads = int(args[4])
    
    data = pd.read_csv(datafile, header=None)
    all_step = data.shape[0] - 1
    
    if not os.path.isdir("result"):
        os.mkdir("result")
    
    for i in tqdm(range(all_step)):
        if threads > 1:
            subprocess.run(['./pci_cluster_dim_ex_parallel.exe', datafile, sigmafile, xifile, str(i), str(threads)])
        else:
            subprocess.run(['./pci_cluster_dim_ex.exe', datafile, sigmafile, xifile, str(i), str(1)])