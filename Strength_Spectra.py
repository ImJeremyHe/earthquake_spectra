#coding=gbk
import os
import os.path
from sdof_nonlinear import *
from read_record import read_record
import os.path

'''
绘制文件夹中的PEER记录的等强度谱
'''


def Strength_Spectra(file_name, t_list,fy_, eta):
    info = read_record(file_name)
    acc = []
    vel = []
    disp = []
    for t in t_list:
        result_dict = Nonlinear_SDOF_Analysis(info, t, 0.05, 10000000000000, eta)
        vme = MASS * result_dict['Ab_Amax']
        fy = fy_ * vme
        result = Nonlinear_SDOF_Analysis(info, t, 0.05, fy, eta)
        acc.append(result['Ab_Amax'])
        vel.append(result['Vmax'])
        disp.append(result['Dmax'])
    return {'t': t_list, 'acc': acc, 'vel': vel, 'disp': disp}
        



if __name__ == '__main__':
    a = Strength_Spectra(r'RSN1_HELENA.A_A-HMC180.AT2',[i/10 for i in range(1, 61)], 0.5, 0)
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
    
