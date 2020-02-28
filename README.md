# Post Hierarchical Clustering Inference
<div align="center">

![](figs/dendrogram_p_black.svg)

</div>

## Abstract
It is an important task to analyze data having multiple clusters behind such as gene expression level data and customer’s purchase history data and discover trends peculiar to each cluster. However, the clustering results are based on subjectivity such as technical knowledge of data, and are not objective. Therefore, we consider using a statistical hypothesis testing to evaluate the reliability of clustering. However, when performing two step inference, such as inference after clustering, the influence of clustering must be considered and corrected appropriately. In this study, we first apply Ward method, which is one of hierarchical clustering, to data with multiple average structures for each cluster to obtain the cluster hierarchical structure. After that, we perform two valid hypothesis testing methods in each branch by exploiting the framework of Selective Inference. 
## Environmental Requirement
- gcc 8.2.0
- GNU Make 3.81
- Install eigen and openmp if compiling c++ source
- Python version 3.7.0
- Please install required packages when the python "ImportError" occurs
  

## Usage
### Hypothesis testing for differences between cluster centers
Under the `cluster` directory
#### (When necessary) Compile 
    
    $ cd /cluster/cpp_source
    $ make
Warning: Change the reference path of `eigen` as needed

Two files are output `pci_cluster_ex`, `pci_cluster_ex_parallel`

#### Preprocessing (When applied to real data)
`preprocess.py`
- Split the data for variance estimation and p-value calculation (variance estimation: p-value calculation = 2: 8)
- Data for p-value calculation is normalized each variables

1. Execute `preprocess.py` for the data you want to preprocess
2. `data`,`stat`, and `interval` directorys are created, and the following files are created in the `data` directory
    - `data.csv` : For p-value calculation
    - `d_ind.csv` : Index which data to use for p-value calculations 
    - `estimate.csv` : For estimation variances
    - `sigma.csv`, `xi.csv` : $\Sigma, \xi$

Example <br>
    
    $ preprocess.py data.csv


#### Clustering and Hypothesis Testing 

##### Calculate p-value for one step
run `calc_p.py` <br>
arguments
- path of `data.csv`
- path of `sigma.csv`
- path of `xi.csv`
- step (0 ~ n - 2)
- Whether to compute in parallel (2 or more: parallel)

Example (parallel computation using 3 cores in the first step) 

    $ pyton calc_p.py data/data.csv data/sigma.csv data/xi.csv 0 3

#### Calculate p-value for all steps
run `calc_p_all.py`
arguments
- path of `data.csv`
- path of `sigma.csv`
- path of `xi.csv`
- Whether to compute in parallel (2 or more: parallel)


Example (not parallel computation) 
 <br>

    $ pyton calc_p_all.py data/data.csv data/sigma.csv data/xi.csv 1


<font color="red"> warning: The file extension of the execute file used in `calc_p_all.py`, `calc_p.py` is `exe`, so please change it appropriately </font>

Both `calc_p.py` and `calc_p_all.py`, the p-value calculation result is output under the `result` directory
- naive_p
- selective_p

`stat` directory is output  a csv file that describes the statistic and the number of dimensions of data at each step, and `interval` directory is output the interval required when calculating the selective-p value <br>

And, `cluster_result` directory is output the following csv file.
- output.csv (It has the same format as `Z` of scikit-learn [linkage](https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.cluster.step.linkage.html))

#### demo
Simple demo data is placed directly under `cluster/data/`

<div align="center">

![50%](figs/demo_data.svg)

</div>

Example <br>

    $ python calc_p_all.py data/demo_data.csv data/demo_sigma.csv data/demo_xi.csv 1

Display dendrogram with p-value <br>

    $ python display_dendro_p_cluster.py

<div align="center">

![50%](figs/demo_dendrogram_p.svg)

</div>


#### demo (FPR)
FPR calculation program in the `demo_FPR` directory <br>
Please run `calc_FPR.py`


Arguments
- number of epoch (default:1000)
- Sample size$n$: ($10, 20, \ldots, n$)
- Dimension $d$
- Whether to compute in parallel (2 or more: parallel)



Example1 ($n = 10, \ldots, 50$, $d = 5$, calculation $1000$ times not parallel)

    $ python calc_FPR.py 1000 50 5 1

例2 (parallel in 4 cores)

    $ python calc_FPR.py 1000 50 5 4


### Hypothesis testing for differences between each dimension of cluster centers
`each_dim`ディレクトリ以下にあるものが相当する

#### Preprocessing (実データ適用時)

- 分散推定用とp値計算用にデータを分ける (分散推定: p値計算 = 2 : 8)
- p値計算用は各変数で正規化される

1. 前処理したいデータに対して, preprocess.pyを実行する
2. `data`のディレクトリが作成され, `data`ディレクトリに次のファイルが作成される
    - `data.csv` : p値計算用
    - `d_ind.csv` : どのデータをp値計算に用いるか 
    - `estimate.csv` : 分散推定用
    - `Sigma.csv`, `Xi.csv` : $\Sigma, \Xi$それぞれ

使用例 <br>
    
    $ preprocess.py data.csv


#### クラスタリングと検定の実行

##### 一つの階層でのp値を計算したい場合
`execute.py`を実行する
以下の引数が必要である
- データファイルパス
- $\Sigma$のファイルパス
- $\Xi$のファイルパス
- どの階層か
- 並列計算するかどうか (2以上で並列となり, 指定したコア数で並列化する)

  

使用例: `$ pyton execute.py data/data.csv data/Sigma.csv data/Xi.csv 0 1` 

#### 全ての階層でのp値を計算したい場合
`execute_allstep.py`を実行する
必要な引数
- データファイルパス
- $\Sigma$のファイルパス
- $\Xi$のファイルパス
- 並列計算するかどうか (2以上で並列となり, 指定したコア数で並列化する)

使用例 <br>

    $ pyton execute_allstep.py data/data.csv data/Sigma.csv data/Xi.csv 1

`execute.py`, `execute_allstep.py`両方とも, 次のcsvファイルが`result`ディレクトリ以下に出力される
- output (クラスタリング結果 scikit-learnのlinkageの出力`Z`と同じフォーマット)
p-valueは各csvファイルに出力される
- naive_p 
- selective_p

#### デモ
クラスタ中心間の差の検定と同様に, `each_dim/data/`直下に簡単なdemoデータが置いてある <br>
1次元目(横軸)ではわかれおらず, 2次元目(縦軸)でわかれているデータ
<div align="center">

![50%](figs/demo_data.svg)

</div>

例 <br>

    $ python execute_allstep.py data/demo_data.csv data/demo_Sigma.csv data/demo_Xi.csv 1

さらに <br>

    $ python display_dendro_p_dim.py 0
と実行すると1次元目における検定でのp値を付与したデンドログラムが得られる. 1次元目ではわかれていないので, 大きな値が得られる.

<div align="center">

![50%](figs/demo_dendrogram_p_dim0.svg)

</div>
縦軸での検定結果. 実際にわかれているので, 一番上の階層においてp値が小さくなる

    $ python display_dendro_p_dim.py 1

<div align="center">

![50%](figs/demo_dendrogram_p_dim1.svg)

</div>


### デンドログラムにp値を付与した図の表示
`/cluster/display_dendro_p_cluster.py`もしくは`/each_dim/display_dendro_p_cluster.py`の`pv_dendrogram`の関数をimportして使ってください <br> <br>
#### pv_dendrogram(sp, nap, start, output, root=0, width=100, height=0, decimal_place=3, font_size=15, **kwargs) 
引数: <br>
- sp: ndarray <br>
    selective-p, ndarrayの次元は1次元にしてください
- nap: ndarray <br>
    naive-p, spと同様に1次元配列にしてください
- start: int <br>
    どの階層からp値を表示するか. 0とすれば, 最初から最後まで表示される. 
- output: list, ndarray <br>
    [scipy.cluster.step.linkage](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.step.linkage.html#scipy.cluster.step.linkage)の`Z`の形式      
- root: int <br>
    デンドログラムの出力を見やすくするために, `Z`の距離に指定された回数rootをとる. 初期値は0.
- width: double <br>
    各階層のnaive-p値とselective-p値を表示する際の横幅. 大きくすればするほど広がる. 初期値100.
- height: double <br>
    各階層のnaive-p値とselective-p値を表示する際の縦幅. 大きくすればするほど上に表示される. 初期値0.
- decimal_place: int <br>
    少数点何桁目まで表示するか. 初期値は3.
- font_size: int <br>
    naive, selectiveのp値および, 凡例の文字サイズ
- **kwargs: <br>
    [scipy.cluster.step.dendrogram](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.step.dendrogram.html#scipy.cluster.step.dendrogram)のkwargsを指定できます.
        
返り値:  <br>
    [scipy.cluster.step.dendrogram](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.step.dendrogram.html#scipy.cluster.step.dendrogram)の出力
    

## ディレクトリ

```
root/
    |- cluster/
          |- cluster_result/
          |- cpp_source/
          |- data/
          |- interval/
          |- stat/
          |- calc_p_all.py
          |- calc_p.py
          |- ...
    |- each_dim/
          |- cpp_source/
          |- data/
          |- result/...
          |- execute_allstep.py
          |- execute.py
          |- ...
    |- figs
    |- README.md
```
### 注意事項
#### dataについて
- 基本的に`data.csv`, `sigma.csv`, `xi.csv`等は値のみのformatとしてください. そうでない場合, errorとなります.

