#coding=gbk
import os
import os.path
from sdof_nonlinear import *
from read_record import read_record
import os.path

'''
绘制文件夹中的PEER记录的等强度谱
'''
def trapezium_area(start_point, end_point, dt):
    return (start_point + end_point) * dt / 2

def Strength_Energy_Spectra(file_name, t_list,fy_, eta, ksai):
    info = read_record(file_name)
    ground_acc = info['ag']
    dt = info['dt']
    energy = []
    for t in t_list:
        result_dict = Nonlinear_SDOF_Analysis(info, t, ksai, 10000000000000, eta)
        vme = 1e-3 * result_dict['Ab_Amax']
        fy = fy_ * vme
        result = Nonlinear_SDOF_Analysis(info, t, ksai, fy, eta)
        e = e_peak = 0
        for i in range(len(ground_acc)-1):
                        e += trapezium_area(-MASS * ground_acc[i] * result['vel'][i],
                                        -MASS * ground_acc[i+1] * result['vel'][i+1], dt)
                        if e > e_peak:
                            e_peak = e
        energy.append(e_peak)
        
    return {'t': t_list, 'energy': energy}
        



if __name__ == '__main__':
    a = Strength_Energy_Spectra(r'RSN1_HELENA.A_A-HMC180.AT2',[i/10 for i in range(1, 61)], 0.5, 0, 0.05)
    with open('energy.txt', 'w') as f:
        for t, energy in zip(a['t'], a['energy']):
            f.write('%-5s%-20s' % (t, energy))
            f.write('\n')

    
