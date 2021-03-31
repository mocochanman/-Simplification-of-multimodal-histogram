import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def weak(y):
    weak_f = []
    for i in range(6):
        weak_f.append(y[i])
    for i in range(6, len(y)):
        ave = y[i]+y[i-1]+y[i-2]+y[i-3]+y[i-4]+y[i-5]+y[i-6]
        weak_f.append(ave/7)    #誤差を減らすため。/で計算
    
    return weak_f

def line_fitted(x, y):

    def gauss_func(x, a, myu, sigma):  #ガウス分布の確率密度関数を定義するお
        return a*x* np.exp( - (np.log(x) - myu)**2 / (2 * sigma**2))

    def line(x, a, b):
        return b + a*x

    y0 = y
    x = np.array([i for i in range(len(y0))])

    #両端点に直線を引いてみよう
    y_end = np.array([y0[0], y0[-1]])
    x_end = np.array([x[0], x[-1]])

    #ここで直線のパラメータを生成
    pop_line, cov = curve_fit(line, x_end, y_end, maxfev=10000000)
    #pop_line=[1/3, 2]
    y_line = []
    for i in x:
        y_line.append(line(i, *pop_line))
    #plt.plot(y0)
    #plt.plot(y_line)

    new_func = []
    for (a,b) in zip(y0, y_line):
        if (a-b) < 0:
            new_func.append(0)
        else:
            new_func.append(a-b)
    new_func = np.array(new_func)
    #plt.plot(new_func)

    popex, cov = curve_fit(gauss_func, x, new_func, maxfev=10000000)
    gauss_new = gauss_func(x, *popex)
    #plt.plot(gauss_new)

    curve_ans = []
    for (a,b) in zip(y_line, gauss_new):
        curve_ans.append(a+b)
    
    return curve_ans