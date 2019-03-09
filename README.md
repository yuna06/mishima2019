# mishima2019


# Single Cell
Single Cell RNA-seq解析では一細胞のみを解析対象とすることで、細胞集団単位の解析では埋もれてしまっていた少数派細胞の発現状態を検出することができます。そのため細胞の発生過程や細胞系譜(cell lineage)を明らかにすることに向いています。
<a href="https://bonohu.github.io/cellfishing.html">坊農秀雅先生のブログ</a>を参考にしました。
CellFishingでは、scRNA-seqで得たデータに適合するデータを検索することができます。

```
# Install Julia
brew cask install julia
# Install CellFishing
julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/bicycle1885/CellFishing.jl.git"))'

https://github.com/bicycle1885/CellFishing.jl.git

# Run test for CellFishing
julia -e 'using Pkg; Pkg.test("CellFishing")'
# git clone CellFishing
git clone https://github.com/bicycle1885/CellFishing.jl
# Install zstd via Homebrew
brew install -v zstd

# Install required packages in julia
julia -e 'using Pkg; Pkg.add("HDF5")'
julia -e 'using Pkg; Pkg.add("DocOpt")'

# Getting data for the search
curl -O http://bimsbstatic.mdc-berlin.de/rajewsky/PSCA/all_sgete_4GU75.loom.gz
gzip -dc all_sgete_4GU75.loom.gz > Plass2018.dge.loom

# Run CellFishing
./bin/cellfishing build Plass2018.dge.loom
./bin/cellfishing search Plass2018.dge.loom.cf Plass2018.dge.loom >neighbors.tsv

```

<a href = "https://hemberg-lab.github.io/scRNA.seq.course/index.html">scRNA.seq.course</a>をやってみました。

### 3 Processing Raw scRNA-seq Data

##### 3.1 FastQC
single-cell RNA-seqデータ解析においても最初はクオリティチェックを行います。今回は遺伝研で貸していただいたMacBook Airで実行を行なったので、ツールのインストールもやっていきます。
```
conda install -c biobuilds fastqc
```
データセットは(Kolodziejczyk et al. 2015)のmESCのデータセットを使います。
シーケンサーはSMART-seq2、ペアエンドです。
(ERR522959_1.fastq, ERR522959_2.fastq)のデータをwgetでダウンロードします。shareというディレクトリを作って、そこにデータをダウンロードします。

```
mkdir share
cd share
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR522/ERR522959/ERR522959_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR522/ERR522959/ERR522959_2.fastq.gz
gunzip ERR522959_1.fastq.gz
gunzip ERR522959_2.fastq.gz
cd ..
```

FastQCを実行していきます。
```
mkdir fastqc_results
fastqc -o fastqc_results Share/ERR522959_1.fastq Share/ERR522959_2.fastq
```

##### 3.2 Trimming Reads
Trim Galoreというソフトウェアを使って不要なリードをトリミングしていきます。
不要なリードとはアダプターやクオリティの低いリードのことです。Trim Galoreのインストールから行なっていきます。
```
conda install -c bioconda trim-galore
```
Trim Galoreを実行していきます。
```
mkdir fastqc_trimmed_results
trim_galore --nextera -o fastqc_trimmed_results Share/ERR522959_1.fastq Share/ERR522959_2.fastq
```
ドキュメントには次にSTARかkallistoでマッピングするようになっていますが、MacBook AirのメモリではSTARを実行できないので、Kallistoを使っていきます。ところで、STARはアライナーである一方で、Kallistoはpseudo-アライナーでアルゴリズムが異なっています。また、STARがリファレンスゲノムにマッピングするのに対し、kallistoではトランスクリプトームにマッピングします。まずはインストールを行います。
```
conda install -c bioconda kallisto
```
##### 3.6.3 Kallisto’s pseudo mode
kallistoのpseudoモードではsingle-cell RNA-seq実験のデータをpseudoアライメントすることができます。
kallitoでは遺伝子ではなくsplice isoformにマッピングしていており、single-cell RNA-seqにおいて次のような理由で難しさが出てきます。
・Single-cell RNA-seqはbulk RNA-seqよりもcoverageが少ない。
・Single-cell RNA-seqプロトコルは3'coverage biasがあるので5'末端のみが異なるisoform同士の場合、判別できない。
・ショートリードのみを使うものではisoform判別が難しくなる。

kallistoを実際に使っていきます。インデックスを作るためのトランスクリプトFASTAは<a href="https://github.com/hemberg-lab/scRNA.seq.course/blob/master/2000_reference.transcripts.fa">こちら</a>からダウンロードできます。
もしくは、hemberg-labの<a href="https://github.com/hemberg-lab/scRNA.seq.course.git">GitHubレポジトリ</a>をクローンしてもいいです。
```
git clone https://github.com/hemberg-lab/scRNA.seq.course.git
```
インデックス作成。
```
mkdir indices
mkdir indices/Kallisto
kallisto index -i indices/Kallisto/transcripts.idx scRNA.seq.course/2000_reference.transcripts.fa
```


```
mkdir results
mkdir results/Kallisto
kallisto pseudo -i indices/Kallisto/transcripts.idx -o results/Kallisto -b batch.txt
```
実行が完了すると、matrix.cells, matrix.ec, matrix.tsv, run_info.jsonというファイルが生成されます。
matrix.ecには一列目にequivalence class IDが二列目にはtranscript IDが書かれています。


# RNA-seq



# 講義
テーブルに酸素濃度の情報などを追加すると、特定のテーマに関するメタアナリシスができる。


# Anacondaアンインストール
Anacondaでこれまで様々なライブラリをインストールしてきたので、ScanpyやCellFishingなど、うまく動かなくなってきました。そこで、今回は一度Anacondaを

```
conda install anaconda-clean
anaconda-clean
cd /Users/username/anaconda3/
rm -rf anaconda3
```

以上でAnacondaのアンインストール、ディレクトリの削除が完了です。
[こちらの記事](https://qiita.com/neet-AI/items/54f675dbcc1a0a513882)に従ってAnacondaをもう一度入れていきます。
```
brew -v
>> Homebrew 2.0.3
brew update
```
公式ホームページからAnacondaダウンロードします。これを機にPython3.6 -> Python3.7にします。
```
conda install python=3.7
```
[condaのチャネルについて](https://qiita.com/yuj/items/8ce25959427ea97d373b)
[jupyter labインストール](https://qiita.com/taka4sato/items/378782763dec3dacb1ee)
[jupyter lab](https://qiita.com/xshirade/items/ca65133aed51d40b5300)
```
conda install -c conda-forge xonsh
conda config --add channels bioconda
# チャンネル一覧を取得
conda config --get channels
conda install -c bioconda scanpy
conda install -c conda-forge jupyterlab
jupyter notebook --generate-config
jupyter lab
```
http://localhost:8080/lab
で表示



# MEMO
- GGGenome
- GGRNA \#
- Metascape
- AOE(All of gene Expression)
- CellFishing
- PCA - 主成分分析
- tSNE \#
- [UMAP](https://umap-learn.readthedocs.io/en/latest/clustering.html) - 次元圧縮
  - [python会記事](https://pythonoum.wordpress.com/2018/09/15/%E8%B6%85%E9%AB%98%E9%80%9F%E6%AC%A1%E5%85%83%E5%9C%A7%E7%B8%AE%E3%82%A2%E3%83%AB%E3%82%B3%E3%82%99%E3%83%AA%E3%82%B9%E3%82%99%E3%83%A0umap/)
- Human Cell Atlas \#
- Tabula Muris \#
- [Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/home) - シングルセルデータを得られる。
- Seurat \#
- [Scanpy](https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html) -
- [anndata](https://anndata.readthedocs.io/en/latest/#)
[python会記事](https://github.com/yyoshiaki/mishima_gassyuku/blob/master/csv2loom/scanpy.ipynb)

\# : added by yyoshiaki

ディレクトリを確認せずに削除
```
rm -rf
```
pigz - マルチコアでgzファイルの圧縮解凍ができる。
HDF - バイナリファイル形式の一つで、一つのファイルの中に階層構造を持たせ、多種類のデータ（文字、数字、行列、画像）を保存することができる。
loompy
