from cProfile import label
from ctypes import BigEndianStructure
from ctypes.wintypes import SIZE
from hashlib import shake_128
from math import log
from re import A, L
from turtle import color
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import re
from scipy.optimize import curve_fit
from typing import Counter
from xml.etree.ElementTree import C14NWriterTarget
from matplotlib import collections
from numpy import int32

rec_con_list=list()
#右側
ref_con_list=list()
#左側
jisin_list=list()
jisbunn2_list=list()
#右in
jisbunn1_list=list()
jisout_list=list()
#左out

rec=list()
ref=list()
space_list=list()


input_file= open("network_judge\\Chemicalreaction_list.txt") 

node_list=list(input_file)
#1つ1つの要素

for u in node_list:
    s=u.split()
    space_list.append(s)
#先ずスペースで分けたリストに変換

for i in range(len(space_list)):
    rec.append(space_list[i][0])
    ref.append(space_list[i][1])
#リストを左右（反応物、生成物）分けて作成

rec_con=Counter(rec)
ref_con=Counter(ref)
#ここの操作
for myvalues in rec_con.values():
    rec_con_list.append(myvalues)


for myvalues2 in ref_con.values():
    ref_con_list.append(myvalues2)


l3=len(rec_con_list)
l4=len(ref_con_list)

c3=Counter(rec_con_list)
c4=Counter(ref_con_list)
#ここで辞書型

for myvalues3 in c3.values():
    jisbunn2_list.append(myvalues3)

for myvalues4 in c4.values():
    jisbunn1_list.append(myvalues4)

for mykey1 in c3.keys():
    jisin_list.append(log(mykey1))


for mykey2 in c4.keys():
    jisout_list.append(log(mykey2))

#リストに戻す


recnode_kind_probab=np.array(jisbunn2_list)#確率
recnode_kind_probab=recnode_kind_probab/l3
recnode_kind_probab=np.log(recnode_kind_probab)
refnode_kind_probab=np.array(jisbunn1_list)
refnode_kind_probab=refnode_kind_probab/l4
refnode_kind_probab=np.log(refnode_kind_probab)

#グラフの縦軸を作成 （各次数）

jisin_list=np.array(jisin_list)#次数計算できるように
jisout_list=np.array(jisout_list)

fig=plt.figure()

ax=fig.add_subplot(2,2,1)#表１k/out用意
ax2=fig.add_subplot(2,2,2)#表2 k/in　用意

coe = np.polyfit(jisin_list,recnode_kind_probab, 1)
print(coe)
coe2=np.polyfit(jisout_list,refnode_kind_probab,1)
print(coe2)

y_fit = coe[0] *jisin_list  + coe[1]
y_fit2 = coe2[0] *jisout_list  + coe2[1]

ax.plot(jisin_list,y_fit, label=str(coe[0]), lw=1,color = "red")
ax2.plot(jisout_list,y_fit2, label=str(coe2[0]), lw=1,color = "red")
#曲線あてはめと、一次式の全体の傾きの出力　coe[0]が傾き
ax.legend()
ax2.legend() #グラフの位置調整

ax.scatter(jisin_list,recnode_kind_probab)
ax2.scatter(jisout_list,refnode_kind_probab)
ax.plot()
ax2.plot()
#中身表示、出力
ax.set_xlabel("k/out")
ax.set_ylabel("p(k/out)")
ax2.set_xlabel("K/in")
ax2.set_ylabel("p(k/in)")
#ラベルを追加
fig.tight_layout()
#詰める
plt.show()

rec_set=set(rec)
ref_set=set(ref)
#セットに変換し、重複を削除
s10_union=rec_set|ref_set
#重複がない状態で結合

#有向グラフのヒストグラム　44 139 ネットワークX
fig2=plt.figure()
ax3=fig2.add_subplot(2,2,1)#表１k/out用意
ax4=fig2.add_subplot(2,2,2)#表2 k/in　用意

ax3.set_yscale("log")
ax3.set_xscale("log")
ax4.set_yscale("log")
ax4.set_xscale("log")
#スケールの変更
G=nx.DiGraph()
#グラフの枠組み
G.add_nodes_from(s10_union)
#中身　化合物
G.add_edges_from(space_list)
#リンク

degree_sequence= sorted(G.in_degree())
degree1_sequence= sorted(G.out_degree())

plt.xscale("log")
plt.yscale("log")

bins =range(0,len(sorted(jisin_list)))
bens=range(0,len(sorted(jisout_list)))

hist,bin=np.histogram(degree_sequence,bins)
hist1,ben=np.histogram(degree1_sequence,bens)
#ヒストグラムを作成
ax3.plot(hist,color="blue")
ax4.plot(hist1,color="red") #赤:out青:in

ax3.set_title('in')
ax4.set_title('out')
plt.show()

input_file.close()











