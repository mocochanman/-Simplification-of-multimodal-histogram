import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import pandas as pd
import sys
import math
import time

#import self definition program
import optimal_pyramid as defi

#l2二乗距離と近似ヒストグラムの計算
def main_algo(right, x, f, prefix_f, ff):
    wl, wu, = defi.weight(right, x, f, prefix_f)
    p = defi.search_peak(wl, wu, right, ff)
    phi_lx, phi_l, phi_ux, phi_u = defi.phi(right, p, prefix_f)
    phi = defi.approximation(right, phi_lx, phi_l, phi_ux, phi_u, f)
    l2 = defi.l2norm(f, phi)
    return l2, phi, p

name = input("file name >>> ")

data = pd.read_csv("/Users/noguchimakoto/study/opt_hist/csv_data/"+str(name)+".csv")
f_origin = list(data['cases'])
f_origin = np.array(f_origin)

#weakly average code
weak_f = []
for i in range(6):
    weak_f.append(f_origin[i])

for i in range(6, len(f_origin)):
    ave = f_origin[i]+f_origin[i-1]+f_origin[i-2]+f_origin[i-3]+f_origin[i-4]+f_origin[i-5]+f_origin[i-6]
    weak_f.append(ave//7)

f_origin = np.array(weak_f)
#end of weakly average

x_plot = np.array([i for i in range(len(f_origin))])
min_id = []


eps = float(input("epsilon >>> "))

for i in range(1):
    f = f_origin
    mt = 0
    count = 0

    find_start = time.time()
    minid = signal.argrelmin(f, order=1)
    find_time = time.time() - find_start

    minimal = minid[0].tolist()
    minimal.insert(0,0)

    appro_phi = []
    left,right,tmp = 0,1,0

    start = time.time()
    prefix = defi.prefix(f)

    while True:
        if len(minimal)==1:
            break
        if right >= len(minimal):    #もし二倍探索が超えたら？
            #まず最後の頂点で試す
            x = [i for i in range(len(f))]
            l2, phi, p = main_algo(len(f)-1, x, f,prefix, f)
            if l2 <= eps:
                appro_phi = appro_phi + phi
                f, prefix = [], []
                break
            else:
                right = len(minimal)-1

        x = [i for i in range(len(f))]
        l2, phi, p = main_algo(minimal[right]-tmp, x, f, prefix, f)

        if right == len(minimal)-1 and l2 <= eps:
            appro_phi = appro_phi + phi
            f, prefix = f[minimal[right]-tmp+1:], prefix[minimal[right]-tmp+1:]
            min_id.append(minimal[right])
            tmp += minimal[left]-tmp+1
            left, right = 0, 1
            mt += 1
            break

        if l2 > eps:    #ここから二分探索
            while(right-left)>1:
                mid = (right+left)//2
                l2, phi, p = main_algo(minimal[mid]-tmp, x, f, prefix, f)
                if eps >= l2:
                   left = mid
                else:
                    right = mid
            l2, phi, p = main_algo(minimal[left]-tmp, x, f, prefix, f)
            appro_phi = appro_phi + phi

            #関数を切って調整する
            min_id.append(minimal[left])
            f, prefix = f[minimal[left]-tmp+1:], prefix[minimal[left]-tmp+1:]
            tmp += minimal[left]-tmp+1
            minimal = minimal[left:]
            left, right = 0, 1
            mt += 1
        else:
            left = right
            right *= 2

    #最後の処理
    if len(f)!=0:
        x = [i for i in range(len(f))]
        l2, phi, p = main_algo(len(f)-1, x, f, prefix, f_origin)
        appro_phi = appro_phi + phi

    elapsed_time = time.time() - start

    l2_2 = defi.l2norm(f_origin, appro_phi)

    print('eps:{} mt:{} L2:{} time:{}'.format(eps, mt+1, l2_2, elapsed_time+find_time))
    plt.bar(x_plot, f_origin,color='blue',edgecolor='black', width=1,linewidth=0.5,label='f')
    #plt.plot(x_plot, f_origin, 'b-', label='f')
    x_phi = [i for i in range(len(appro_phi))]
    appro_phi = np.array(appro_phi)
    plt.title('New infections of covid-19 in Japan')
    plt.xlabel('x')
    plt.ylabel('y')
    #plt.plot(x_phi, appro_phi, 'r-', label='φ')
    #plt.bar(x_phi, appro_phi, color='yellow',edgecolor='black', width=1,linewidth=1, alpha=0.5,label='φ')
    #plt.scatter(x_plot[min_id],appro_phi[min_id], color='red', zorder=3,label='deviding point')
    plt.legend()
    plt.show()