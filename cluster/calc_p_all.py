import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import sys
import subprocess
import os
from tqdm import tqdm
from mpmath import *


def calc_norm_mp(stat, final_interval, digit=1000):
    """
    高精度な切断正規分布でのp値計算
    """
    mp.dps = digit
    l = final_interval[:, 0]
    u = final_interval[:, 1]
    deno = 0
    nume = 0
    for i in range(len(l)):
        deno += mp.ncdf(u[i]) - mp.ncdf(l[i])
        if l[i] < stat < u[i]:
            nume += mp.ncdf(u[i]) - mp.ncdf(stat)
        elif u[i] > stat :
            nume += mp.ncdf(u[i]) - mp.ncdf(l[i])
    return float(nume / deno)

def chi2cdf(stat, df, digit=1000):
    """
    カイ二乗分布の累積分布関数
    """
    mp.dps = digit
    if stat < 0:
        return 0.0
    else:
        gamma = mp.gamma(float(df) / 2)
        imcomp = mp.gammainc(float(df) / 2, 0, float(stat) / 2)
        return imcomp / gamma
    
def calc_chi2_mp(stat, df, final_interval, digit=1000):
    """
    高精度な切断カイ二乗分布でのp値計算
    """
    l = final_interval[:, 0]
    u = final_interval[:, 1]
    deno = 0
    nume = 0
    for i in range(len(l)):
        deno += chi2cdf(u[i], df, digit) - chi2cdf(l[i], df, digit)
        if l[i] < chi2 < u[i]:
            nume += chi2cdf(u[i], df, digit) - chi2cdf(chi2, df, digit)
        elif u[i] > chi2 :
            nume += chi2cdf(u[i], df, digit) - chi2cdf(l[i], df, digit)
    return float(nume / deno)

if __name__ == '__main__':
    
    args = sys.argv
    datafile = args[1]
    sigmafile = args[2]
    xifile = args[3]
    all_step = int(args[4])
    
    xi = pd.read_csv(xifile, header=None).values.reshape(-1, )[0]
    
    if not os.path.isdir("cluster_result"):
        os.mkdir("cluster_result")
    
    naive_p_all = []
    selective_p_all = []
    selective_p_aprox_all = []
    for i in tqdm(range(all_step)):
        statfile = "stat/test_stat" + "_step" + str(i) + ".csv"
        intervalfile = "interval/final_interval" + "_step" + str(i) + ".csv"
        if os.path.exists(statfile):
            os.remove(statfile)
        if os.path.exists(intervalfile):
            os.remove(intervalfile)
            
        subprocess.run(['./pci_cluster_ex.exe', datafile, sigmafile, str(xi), str(i)])

        # 結果 読み込み
        chi2, d = pd.read_csv(statfile, header=None).values.reshape(2, )
        final_interval = pd.read_csv(intervalfile, header=None).values.reshape(-1, 2)

        naive_p = 1 - stats.chi2.cdf(chi2, d)
        selective_p = calc_chi2_mp(chi2, d, final_interval)
        
        naive_p_all.append(naive_p)
        selective_p_all.append(selective_p)

    naive_p_all = pd.DataFrame(naive_p_all)
    selective_p_all = pd.DataFrame(selective_p_all)
    
    if not os.path.isdir("result"):
        os.mkdir("result")
    naive_p_all.to_csv("result/naive_p.csv", header=None, index=False)
    selective_p_all.to_csv("result/selective_p.csv", header=None, index=False)