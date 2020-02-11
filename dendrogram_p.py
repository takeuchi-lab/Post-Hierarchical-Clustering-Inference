import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, fcluster

def pv_dendrogram(sp, nap, start, output, root=0, width=100, height=0, decimal_place=3, font_size=15, **kwargs):
    """
    display dendrogram with selective-p and naive-p
    args
    ----------
    sp: ndarray
        selective-p
    nap: ndarray
        naive-p
    start: int
        From which step to display
    output: list, ndarray
        Z that is the output of scipy.cluster.hierarchy.linkage: list or ndarray
    kwargs
    ----------
    root: int
        To make the output of the figure easier to see, the root of the distance is taken root times
    width: double
        how far apart naive-p and selective-p
    height: double
        height to display naive-p and selective-p
    decimal_place: int
        Number of decimal places
    font_size: int
        fontsize of naive-p, selective-p and legend
    **kwargs: 
        **kwargs of scipy.cluster.hierarchy.dendrogram
        
    Returns
    -------
        output of scipy.cluster.hierarchy.dendrogram
    """    
    
    stepsize = output.shape[0]
    
    output = np.array(output)
    for i in range(root):
        output[:, 2] = np.sqrt(output[:, 2])
    
    ddata = dendrogram(output, **kwargs)
    
    xarray = np.array(ddata['icoord'])
    yarray = np.array(ddata['dcoord'])
    
    xarray = xarray[np.argsort(yarray, axis=0)[:, 2]]
    yarray = yarray[np.argsort(yarray, axis=0)[:, 2]]
    
    for i in range(stepsize):
        if i >= start:
            xm = 0.5 * (xarray[i][1] + xarray[i][2])
            xdiv = abs(xarray[i][1] + xarray[i][2])
            # x1, x2の幅は文字の大きさや図の大きさによって適切に変更する必要がある
            x1 = xm - width * (xm/xdiv)
            x2 = xm + width * (xm/xdiv)

            y1 = yarray[i][1] + height
            y2 = yarray[i][1] + height
            
            plt.annotate("{}".format(np.round(sp[i], decimal_place)), (x2, y2), color="fuchsia", xytext=(0, -5), textcoords='offset points', va='top', ha='center', label="s", fontsize=font_size)
            plt.annotate("{}".format(np.round(nap[i], decimal_place)), (x1, y1), color="darkcyan", xytext=(0, -5), textcoords='offset points', va='top', ha='center', label="n", fontsize=font_size)
        
    plt.annotate("selective", color='magenta',xy=(np.max(xarray) - np.max(xarray)*0.15, max(output[:, 2])), fontsize=font_size)
    plt.annotate("naive", color='darkcyan',xy=(np.max(xarray) - np.max(xarray)*0.2, max(output[:, 2])), fontsize=font_size)
    return ddata

"""
example
"""
if __name__ == "__main__":
    
    step = 5
    
    output = pd.read_csv("result/output.csv", header=None).iloc[:, :4].values
    
    sp = pd.read_csv("result/selective-p.csv", header=None).values
    nap = pd.read_csv("result/naive-p.csv", header=None).values
    
    pv_dendrogram(sp, nap, 0, output, root=0, width=20, height=0, decimal_place=3, font_size=15, color_threshold=5, leaf_font_size=15)
    plt.savefig("dendrogram_p.svg", pad_inches=0, bbox_inches='tight')
    plt.show()