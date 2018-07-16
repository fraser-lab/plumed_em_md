>> EMDB, PDB 用のGMMの作成・更新方針<<
2015年7月15日(水)
川端　猛

========EMDB用のGMMの更新========

>> EMDB のデータの更新　<<
(1) dodeca にログイン
(2) cd DB/emdb/structures
(3) ./pdbj_update　　　　※このアップデートが結構時間がかかる。３０分ぐらい？
>> GMMの計算 <<
(4) cd /home/takawaba/etc/EMDB_GMM_CURRENT/
(5) rm UPDATE_EMDB_GMM_NG20/*
(6) ~/work/gmconvert/src/emdb_python/emdb2gmm.py -ogdir EMDB_GMM_NG20 -ng 20 -update T -ogdirup UPDATE_EMDB_GMM_NG20 -A T
    ※更新数が多い場合は、必要に応じて、-div [bunshi]/[bunbo] オプションをつける。
>> GMMの転送　<<
(7) ssh kawabata@pdbjkw1.pdbj.org
(8) cd /var/www/html/gmfit
(9) cd EMDB_GMM_NG20
(10) sftp takawaba@192.168.39.136
(11) sftp> cd etc/EMDB_GMM_CURRENT/UPDATE_EMDB_GMM_NG20/
(12) sftp> get *.gmm
(13) sftp> quit

======= PDB用のGMMの更新===========
>> mmCIFのPDBデータの更新 <<
(1) dodecaにログイン
(2) cd DB/mmCIF
(3) ./pdbj_update
>> GMMの計算 <<
(4) cd /home/takawaba/etc/PDB_GMM_CURRENT/
(5) rm -r UPDATE_PDB_GMM_NG20/PDB_GMM_NG20/*
(6) ~/work/gmconvert/src/emdb_python/mmCIF2gmm.py -ng 20 -odir PDB_GMM_NG20 -update T -odirup UPDATE_PDB_GMM_NG20 -A 
 ※更新数が多い場合は、必要に応じて、-div [bunshi]/[bunbo] オプションをつける。
(7) cd UPDATE_PDB_GMM_NG20
(8)  tar cvf UPDATE_PDB_GMM_NG20_[date].tar PDB_GMM_NG20

>> GMMの転送　<<
(9) ssh kawabata@pdbjkw1.pdbj.org
(10) cd /var/www/html/gmfit
(11) sftp takawaba@192.168.39.136
(12) sftp> cd etc/PDB_GMM_CURRENT/UPDATE_PDB_GMM_NG20
(13) sftp> get UPDATE_PDB_GMM_NG20_[date].tar
(14) sftp> quit
(15) tar xvf UPDATE_PDB_GMM_NG20_[date].tar




[EMDB用のGMMの作成・更新について]

EMDBの密度マップを表現するのに必要なgdfの数は、
格子数と解像度の両方に依存するためなかなか難しい。

よって、現状ではすべて一律の数のgdfを設定することにする。
どのくらいの数が適当であるかは、いろいろな考え方によるが、
例えばウィルスなど正20面体の対称性を持つ場合、少なくとも20は必要だろう。
正20面体の１面が三つのサブユニットで構成されていた場合は60が必要ということになる。

とりあえず、現状では40のgdfを一律に用いることにしようと、2015/06/20では
思っていたが、omokageの投稿論文では、一律Ngauss=20と書いてしまったので、
しばらく、Ngauss=20で行うことにする。


GMMを作成する計算はdodecaで行う。

/home/takawaba/etc/EMDB_GMM_CURRENTというディレクトリの下の、
EMDB_GMM_NG40に作成したGMMを収納する。
ここで、毎回の更新時のときの新しく作成した分だけ、
UPDATE_EMDB_GMM_NG40というディレクトリにコピーし、
これをtarでまとめて、サーバマシンpdbjkw1.pdbj.orgにコピーするという方針をとる。
例えば、以下のコマンドで作成する。

~/work/gmconvert/src/emdb_python/emdb2gmm.py -ogdir EMDB_GMM_NG20 -ng 20 -update T -ogdirup UPDATE_EMDB_GMM_NG20 -A T 

~/work/gmconvert/src/emdb_python/emdb2gmm.py -ogdir EMDB_GMM_NG40 -ng 40 -update T -ogdirup UPDATE_EMDB_GMM_NG40 -A T 

もし、数が多いようなら　-div 0/4,... などのオプションをつけて、
複数個のジョブを並列で実行すればよい。

[PDB用のGMMの作成・更新について]
現行では、PDBについては、二重の構造になっている。一つは
(1)PDBのasymmetric unitとbiological unitを一定のNgdf=20 or 40で表現する
(2)PDBの各鎖を残基数に応じた数のgdfで表現する方法
の二つである。

残基に応じて、Ngaussを決めるには、１ドメイン５０アミノ酸を最低２つのgdfで表現することを
考えると、

Ngauss = ceil[(Naa-50)/25] + 2 

という関数になる。この算法によれば、Naa=500で、Ngauss=20, Naa=1000で、Ngauss=40となる。
ダイニン(3vkhA)は3042アミノ酸、リボゾーム(4v69)は8000アミノ酸を越える。ウィルスはもっと大きくなる。
よって、1000アミノ酸を越えたら、一律　Ngauss=40というのが妥当な計算になるような気がする。

ウィルスなどが非常に計算時間がかかるので、現状では、PDBについては、一律Ngauss=20ということで
設定してある。

この場合のコマンドは以下のとおり。
nohup ~/work/gmconvert/src/emdb_python/mmCIF2gmm.py -ng 20 -odir PDB_GMM_NG20 -update T -odirup UPDATE_PDB_GMM_NG20 -A T -div 0/8 > 0.log &





