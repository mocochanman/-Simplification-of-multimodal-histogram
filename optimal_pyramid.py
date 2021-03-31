import numpy as np
import math
import matplotlib.pyplot as plt

def prefix(f):
    #fの積分(プレフィックス和)Fを計算する
    prefix = [f[0]]
    for i in range(1, len(f)):
        prefix.append(f[i]+prefix[i-1])

    return prefix

#上下凸包の各頂点までの重みを計算する
def weight(n, x, y, prefix):
    n = n+1
    #凸包木の計算
    l = []      #下側凸包のx座標
    Tl = []     #下側凸包のy座標
    wl = [0]    #各頂点までの凸包の辺の重み

    #下側凸包木
    #初めの２頂点のみ格納
    for i in range(2):
        l.append(x[i])
        Tl.append(prefix[i])
    wl.append((x[1]-x[0])*((prefix[1]-prefix[0])/(x[1]-x[0])) * ((prefix[1]-prefix[0])/(x[1]-x[0])))    #((xi - xi-1)t(e)^2)

    #凸包の計算
    for i in range(2, n):
        #ベクトルの外積を計算し、非凸包であればその頂点を省く（２番目の頂点）
        while len(l)>=2 and ((l[-1]-l[-2])*(prefix[i]-Tl[-2]) - (i-l[-2])*(Tl[-1]-Tl[-2]))<0:
            l.pop()
            Tl.pop()
        #凸になる時の値を各リストに格納していく
        wl.append(wl[l[-1]] + (i-l[-1]) * ((prefix[i]-prefix[l[-1]])/(i-l[-1])) * ((prefix[i]-prefix[l[-1]])/(i-l[-1])))
        l.append(i)
        Tl.append(prefix[i])

    u = []    #上側凸包のx座標
    Tu = []   #上側凸包のy座標
    wu = [0]

    #上側凸包木
    #初めの２頂点のみ格納
    u.append(x[n-1])
    u.append(x[n-2])
    Tu.append(prefix[n-1])
    Tu.append(prefix[n-2])
    wu.append( (x[n-1]-x[n-2]) * ((prefix[n-1]-prefix[n-2])/(x[n-1]-x[n-2]))*((prefix[n-1]-prefix[n-2])/(x[n-1]-x[n-2])))

    #凸包の計算（外積）
    for i in reversed(range(0, n-2)):
        while len(u)>=2 and ((u[-1]-i)*(Tu[-2]-prefix[i]) - (u[-2]-i)*(Tu[-1]-prefix[i]))>=0:
            u.pop()
            Tu.pop()
        wu.append(wu[n-1-u[-1]] + (u[-1]-i) * ((prefix[u[-1]]-prefix[i])/(u[-1]-i)) * ((prefix[u[-1]]-prefix[i])/(u[-1]-i)))
        u.append(i)
        Tu.append(prefix[i])
    wu.reverse()

    return wl, wu


#ピーク頂点を探索する関数
def search_peak(wl, wu, n, ff):
    peak = 0
    max_weight = wl[peak] + wu[peak]    #φの重みの初期値(ピークを０としている)
    for i in range(1, n+1):
        tmp = wl[i] + wu[i]
        if (max_weight < tmp):
            max_weight = tmp
            peak = i
    return peak


#ピーク決定後に近似曲線を計算するプログラム
def phi(n, p, prefix):
    n =  n+1
    #左側の初期化
    l = [0, 1]
    Tl = [prefix[0], prefix[1]]

    #左側の下側凸包を計算
    for i in range(2, p+1):
        while len(l)>=2 and ((l[-1]-l[-2])*(prefix[i]-Tl[-2]) - (i-l[-2])*(Tl[-1]-Tl[-2]))<0:
            l.pop()
            Tl.pop()
        l.append(i)
        Tl.append(prefix[i])

    #右側の初期化
    u = [n-1, n-2]
    Tu = [prefix[n-1],prefix[n-2]]
    for i in reversed(range(p, n-1)):
        while len(u)>=2 and ((u[-1]-i)*(Tu[-2]-prefix[i]) - (u[-2]-i)*(Tu[-1]-prefix[i]))>=0:
            u.pop()
            Tu.pop()
        u.append(i)
        Tu.append(prefix[i])
    
    #ピークインデックスが０の時の処理
    if p == 0:
        return [], [], u, Tu
    elif p == n-1:
        return l, Tl, [], []
    else:
        return l, Tl, u, Tu


#簡略化したヒストグラムを計算するプログラム
def approximation(n, l, Tl, u, Tu, f):
    phi = []
    for i in range(len(l)-1):
        if i == 0:
            if f[1] <= f[0] or l[i+1]-l[i]!=1:
                phi.append((Tl[i+1]-Tl[i])/(l[i+1]-l[i]))
            else:
                phi.append(f[0])
        for j in range(l[i+1]-l[i]):
            phi.append((Tl[i+1]-Tl[i])/(l[i+1]-l[i]))

    if len(l) == 0:
        phi.append(f[0])
    for i in reversed(range(1, len(u))):
        for j in range(u[i-1]-u[i]):
            phi.append((Tu[i-1]-Tu[i])/(u[i-1]-u[i]))

    return phi


#l2二乗距離を計算するプログラム
def l2norm(f, phi):
    l2 = 0
    for i in range(len(phi)):
        l2 += (f[i]-phi[i]) * (f[i]-phi[i])
    return l2
