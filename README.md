# ネットワークカテゴリ判定プログラム（次数分布を作成）

生成物と反応物を示した反応物と生成物のテキストデータから、次数分布を作成してグラフの概要からネットワークカテゴリを判定します。

## プログラムのプロセス

* pythonで各ノードの次数の取得、計算を行う。

* Matplotlib次数分布（散布図）の作成（スケールをlogでとる）

* 曲線あてはめの適用(Curve fitting)

* この図の概形からスケールフリーネットワークかどうか判定する。

* 追加で、networkXによるヒストグラムの頂点を結んだグラフも作成した。