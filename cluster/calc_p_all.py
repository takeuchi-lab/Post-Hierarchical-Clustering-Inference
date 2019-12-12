import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import sys
import subprocess
import os
from tqdm import tqdm

def calc_p(chi2, df, final_interval):
    """
    近似計算なしの切断カイ二乗分布でのp値計算
    """
    l = final_interval[:, 0]
    u = final_interval[:, 1]
    deno = np.sum(stats.chi2.cdf(u,df) - stats.chi2.cdf(l,df))
    p = np.argmax( (l < chi2) * (chi2 < u) )
    add = np.sum( stats.chi2.cdf(u[p+1:], df) - stats.chi2.cdf(l[p+1:], df) ) if p + 1 < len(l) else 0.0
    nume = stats.chi2.cdf(u[p], df) - stats.chi2.cdf(chi2, df) + add 
    return nume / deno

def ImpSamp(sig2, tau, L, U, num_samp=10**6):
    x = np.random.normal(tau, sig2, num_samp)
    arg = (-tau * x / sig2) + (tau ** 2 / sig2)
    ex = np.exp( arg )
    nume = np.sum( ex[(tau <= x) * (x <= U)] )
    deno = np.sum( ex[(L <= x) * (x <= U)] )
    return nume / deno

"""区間は昇順にソートされていることを想定"""
def calc_multi_p(tau, mean, std, L, U):
    deno = np.sum(stats.norm.cdf(-(L - mean) / std) - stats.norm.cdf(-(U - mean) / std))
    if deno <= 0:
        if len(L) < 2:
            return ImpSamp(std, tau, L, U)
        else:
            p1 = ImpSbamp(std, U[0], L[0], U[1])
            p2 = ImpSamp(std, L[1], L[0], U[1])
            p = np.argmax(((L < tau) * (tau < U)))
            if p == 0:
                deno = ( 1 + p2 / (1 - p1) ) + 1. / ( 1 + (1 - p1) / p2 )
                nume = ImpSamp(std, tau, L[0], U[0])
                return nume / deno
            else:
                deno = 1 + (1 - p1) / p2
                nume = ImpSamp(std, tau, L[1], U[1])
                return nume / deno
    else:
        p = np.argmax( (L < tau) * (tau < U) )
        add = np.sum( stats.norm.cdf(-(L[p+1:] - mean) / std) - stats.norm.cdf(-(U[p+1:] - mean) / std) ) if p + 1 < len(L) else 0.0
        nume = stats.norm.cdf(-(tau - mean) / std) - stats.norm.cdf(-(U[p] - mean) / std) + add
        return nume / deno


def calc_p_approx(chi2, df, final_interval):
    """
    切断カイ二乗分布でのp値計算 正規近似var
    """
    L = final_interval[:, 0]
    U = final_interval[:, 1]
    # 近似
    Lx = (chi2 / df)**(1. / 6) - (chi2 / df)**(1. / 3) / 2 + (chi2 / df)**(1. / 2) / 3
    mean = 5. / 6 - 1. / (9 * df) - 7. / (648 * df**2) + 25. / (2187 * df**3)
    var = 1. / (18 * df) + 1. / (162 * df**2) - 37. / (11664 * df**3)
    
    LL = (L / df)**(1. / 6) - (L / df)**(1. / 3) / 2 + (L / df)**(1. / 2) / 3
    if (U != np.inf).all():
        LU = (U / df)**(1. / 6) - (U / df)**(1. / 3) / 2 + (U / df)**(1. / 2) / 3
    elif (U == np.inf).all():
        LU = U
    selective_p = calc_multi_p(Lx, mean, np.sqrt(var), LL, LU)
    return selective_p



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
        selective_p = calc_p(chi2, d, final_interval)
        selective_p_aprox = calc_p_approx(chi2, d, final_interval)
        
        naive_p_all.append(naive_p)
        selective_p_all.append(selective_p)
        selective_p_aprox_all.append(selective_p_aprox)

    naive_p_all = pd.DataFrame(naive_p_all)
    selective_p_all = pd.DataFrame(selective_p_all)
    selective_p_aprox_all = pd.DataFrame(selective_p_aprox_all)
    
    if not os.path.isdir("result"):
        os.mkdir("result")
    naive_p_all.to_csv("result/naive_p.csv", header=None, index=False)
    selective_p_all.to_csv("result/selective_p.csv", header=None, index=False)
    selective_p_aprox_all.to_csv("result/selective_p_aprox.csv", header=None, index=False)