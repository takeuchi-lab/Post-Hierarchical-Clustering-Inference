# Post Hierarchical Clustering Inference
<div align="center">

![](figs/dendrogram_p_black.svg)

</div>

## Abstract
DNA マイクロアレイにより計測された遺伝子発現量データの解析が遺伝学分野で盛んにおこなわれている. 例えば, マイクロアレイデータに対してクラスタリングを適用し, 得られたクラスタを特徴づける遺伝子の統計的有意性を検証する. これにより, 疾患のサブタイプを同定するという試みである. しかし, このようなクラスタリング後の推論において古典的な検定手法を用いてしまうと, クラスタリングの影響を考慮できず, 検定におけるtype I error を制御できない. そこで, 本研究では, 階層型クラスタリングのウォード法に対して, Lee et al. [1] によって提案されたSelective Inference の枠組みを適用し, 妥当な検定を行う手法を提案する.
## Environmental Requirement
- c++ソースをコンパイルする場合はeigenをインストールする
- Python version 3.7.0
- Please install required packages when the python "ImportError" occurs
  

## Usage
### クラスタ中心間の検定
`cluster`ディレクトリ以下にあるものが相当する
#### (必要時) コンパイル 
次のコマンドで実行ファイルが`/cluster`直下に作成される (注意: 必要に応じて, eigenの参照パスを変更する必要がある)
    
    $ cd /cluster/cpp_source
    $ make
`pci_cluster_ex`, および並列実行用`pci_cluster_ex_parallel` の２つ

#### Preprocessing (実データ適用時)

- 分散推定用とp値計算用にデータを分ける (分散推定: p値計算 = 2 : 8)
- p値計算用は各変数で正規化される

1. 前処理したいデータに対して, preprocess.pyを実行する
2. `data`, `stat`, および`interval`のディレクトリが作成され, `data`ディレクトリに次のファイルが作成される
    - `data.csv` : p値計算用
    - `d_ind.csv` : どのデータをp値計算に用いるか 
    - `estimate.csv` : 分散推定用
    - `sigma.csv`, `xi.csv` : $\Sigma, \xi$それぞれ

使用例 <br>
    
    $ preprocess.py data.csv


#### クラスタリングと検定の実行

##### ある階層でのp値を計算したい場合
`calc_p.py`を実行する
以下の引数が必要である
- データファイルパス
- $\Sigma$のファイルパス
- $\xi^2$のファイルパス
- どの階層か (0 ~ n - 2)
- 並列計算するかどうか (2以上で並列となり, 指定したコア数で並列化する)

使用例: <br>

    $ pyton calc_p.py data/data.csv data/sigma.csv data/xi.csv 0 1

#### 全ての階層でのp値を計算したい場合
`calc_p_all.py`を実行する
必要な引数
- データファイルパス
- $\Sigma$のファイルパス
- $\xi^2$のファイルパス
- 並列計算するかどうか (2以上で並列となり, 指定したコア数で並列化する)

使用例 <br>

    $ pyton calc_p_all.py data/data.csv data/sigma.csv data/xi.csv 1


<font color="red"> 注: `calc_p_all.py`, `calc_p.py`で使われる実行ファイルの拡張子は`exe`なので適宜変更してください </font>

`calc_p.py`, `calc_p_all.py`両方とも, `result`ディレクトリ以下にp値計算結果が出力される
- naive_p
- selective_p
  
`stat`には, 各階層での統計量とデータの次元数を記載したcsvファイルが, `interval`にはselective-p値計算時に必要な区間が出力される

また, `cluster_result`ディレクトリ以下に次のものが出力される
- output (クラスタリング結果 scikit-learnの[linkage](https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.cluster.hierarchy.linkage.html)の出力`Z`と同じフォーマット)

#### デモ
`cluster/data/`直下に簡単なdemoデータが置いてある

<div align="center">

![50%](figs/demo_data.svg)

</div>

例 <br>

    $ python calc_p_all.py data/demo_data.csv data/demo_sigma.csv data/demo_xi.csv 1

さらに <br>

    $ python display_dendro_p_cluster.py
と実行するとp値を付与したデンドログラムが得られる

<div align="center">

![50%](figs/demo_dendrogram_p.svg)

</div>


#### デモ (FPR)
`demo_FPR`ディレクトリにFPRの計算プログラムがある <br>
`calc_FPR.py`を実行することで, 計算される



必要な引数
- 繰り返し回数 (ここでは1000回)
- $n$の最大サイズ: ($10, 20, \ldots, n$)
- 次元 $d$
- 並列オプション : 2以上であれば, 並列計算を指定したコア数


例1 ($n = 10, \ldots, 50$, $d = 5$, $1000$回の計算 並列化しない)

    $ python calc_FPR.py 1000 50 5 1

例2 (4コアで並列)

    $ python calc_FPR.py 1000 50 5 4


### クラスタ中心の各次元での検定
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
    [scipy.cluster.hierarchy.linkage](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage)の`Z`の形式      
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
    [scipy.cluster.hierarchy.dendrogram](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.dendrogram.html#scipy.cluster.hierarchy.dendrogram)のkwargsを指定できます.
        
返り値:  <br>
    [scipy.cluster.hierarchy.dendrogram](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.dendrogram.html#scipy.cluster.hierarchy.dendrogram)の出力
    

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
- data.csvは値のみのformatとしてください. そうでない場合, errorとなります.

