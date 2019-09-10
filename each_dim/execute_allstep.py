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
    all_step = int(args[4])
    
    if not os.path.isdir("result"):
        os.mkdir("result")
    
    naive_p_all = []
    selective_p_all = []
    selective_p_aprox_all = []
    for i in tqdm(range(all_step)):
        subprocess.run(['./pci_cluster_step_ex.exe', datafile, sigmafile, xifile, str(i)])