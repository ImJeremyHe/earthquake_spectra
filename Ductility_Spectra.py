#coding=gbk
import os
import os.path
from sdof_nonlinear import *
from read_record import read_record
from math import log, exp

'''
绘制文件夹中的PEER记录的等延性谱
'''


def Ductility_Spectra(file_name, t_list,d, eta, damping_ratio):
    info = read_record(file_name)
    initial = [0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 0.9]
    recorder = []
    acc = []
    vel = []
    disp = []
    for t in t_list:
        result_dict = Nonlinear_SDOF_Analysis(info, t, damping_ratio, 10000000000000, eta)
        vme = 1e-3 * result_dict['Ab_Amax']
        result = iteration(info, t, damping_ratio, eta, d, vme, initial, 3000)
        acc.append(result['Ab_Amax'])
        vel.append(result['Vmax'])
        disp.append(result['Dmax'])
    return {'t': t_list, 'acc': acc, 'vel': vel, 'disp': disp}
        

def iteration(info, t, damping_ratio, eta, d, vme, initial, max_times = 1000):
    ''' recorder - 记录标准屈服强度与延性系数对应关系。(fy_, mu)
        info - 地震动记录信息
        t - 周期
        damping_ratio - 阻尼比
        eta - 第二刚度折减系数
        d - 目标延性

    '''
    recorder = []
    #print("============================================================")
    for fy_ in initial:
        fy = vme * fy_
        result = Nonlinear_SDOF_Analysis(info, t, damping_ratio, fy, eta)
        if abs((result['mu'] - d )/ d) <= 0.05:
            return result
        recorder.append((fy_, result['mu']))
    for i in range(max_times):
        #print('-------------------------------------------------------')
        recorder = sorted(recorder, key=lambda x: abs(x[1]-d))
        a = recorder[0]
        b = recorder[1]
        #print(recorder)
        predict_fy_ = exp((log(a[0]) - log(b[0])) / (log(a[1]) - log(b[1])) * (log(d) - log(a[1])) + log(a[0]))
        fy = vme * predict_fy_
        result = Nonlinear_SDOF_Analysis(info, t, damping_ratio, fy, eta)
        if abs((result['mu'] - d) / d) <= 0.05:
            return result
        del recorder[i%(i%3+1)]
        #print(result['mu'])
        #print(recorder)
        recorder.append((predict_fy_, result['mu']))
    raise Exception


if __name__ == '__main__':
    a = Ductility_Spectra(r'RSN1_HELENA.A_A-HMC180.AT2',[i/10 for i in range(1, 51)], 4, 0, 0.05)
    with open('acc_spect.txt', 'w') as f:
        for t, acc in zip(a['t'], a['acc']):
            f.write('%-5s%-20s' % (t, acc))
            f.write('\n')
    with open('vel_spect.txt', 'w') as f:
        for t, vel in zip(a['t'], a['vel']):
            f.write('%-5s%-20s' %(t, vel))
            f.write('\n')
    with open('disp_spect.txt', 'w') as f:
        for t, disp in zip(a['t'], a['disp']):
            f.write('%-5s%-20s' %(t, disp))
            f.write('\n')
    
