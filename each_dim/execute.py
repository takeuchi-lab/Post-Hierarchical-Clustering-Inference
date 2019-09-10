import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import sys
import subprocess
import os

if __name__ == '__main__':
    
    args = sys.argv
    datafile = args[1]
    sigmafile = args[2]
    xifile = args[3]
    step = args[4]
        
    if not os.path.isdir("result"):
        os.mkdir("result")

    subprocess.run(['./pci_cluster_step_ex.exe', datafile, sigmafile, xifile, step])