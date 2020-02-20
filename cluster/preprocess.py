import numpy as np 
import pandas as pd
import sys
import os

if __name__ == "__main__":
    if not os.path.isdir("data"):
        os.mkdir("data")
    if not os.path.isdir("stat"):
        os.mkdir("stat")
    if not os.path.isdir("interval"):
        os.mkdir("interval")

    args = sys.argv
    datafile = args[1]
    data = pd.read_csv(datafile, header=None)

    n = data.shape[0]
    ind = np.arange(0, n)
    d_ind = np.sort(np.random.choice(ind, int(n*0.8), replace=False))
    e_ind = np.sort(np.array(list(set(ind) - set(d_ind))))

    data_p = data.iloc[d_ind]
    esti = data.iloc[e_ind]

    data_p_preprocess = (data_p - data_p.mean()) / data_p.std()
    esti_preprocess = (esti - esti.mean()) / esti.std()

    data_p_preprocess.to_csv("data/data.csv", header=None, index=False)
    esti_preprocess.to_csv("data/estimate.csv", header=None, index=False)
    pd.DataFrame(d_ind).to_csv("data/d_ind.csv", header=None, index=False)

    # 推定用のデータサイズが実験用のデータサイズと異なる場合
    sigma = np.mean(esti_preprocess.T.var())
    Sigma = sigma * np.identity(data_p_preprocess.shape[0])
    
    Sigma = pd.DataFrame(Sigma)
    Sigma.to_csv("data/sigma.csv", header=None, index=False)
    
    xi = np.mean(esti_preprocess.var())
    xi = pd.DataFrame(np.array([xi]))
    xi.to_csv("data/xi.csv", header=None, index=False)
