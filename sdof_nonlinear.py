from read_record import *
from force_deformation import *
MASS = 1e-3

def Nonlinear_SDOF_Analysis(record_info, T, ksai, strength, eta):
    m = MASS
    k = m / ((T/2/3.14159) ** 2)
    c = ksai * 2.0 * (m*k) ** 0.5
    props = [k, strength, eta]
    Amax = 0
    Vmax = 0
    Dmax = 0
    kt = 0
    ke = 0
    du = 0
    u_trial = 0
    du_trial = 0
    force_c = 0
    force_p = 0
    Res = 0
    df = 0
    statev = [0, 0, 0, 0, 0, 0, 0]

    dt = record_info['dt']
    ag = record_info['ag']

    beta = 0.25
    gama = 0.5
    acc = []
    vel = []
    dis = []
    abs_acc = []
    force = []
    force.append(0.0)
    dis.append(0.0)
    vel.append(0.0)
    acc.append((-ag[0]*m - c*vel[0] - k*dis[0])/m)

    a = m / beta / dt + gama * c / beta
    b = m/2/beta + dt*(gama / 2.0 / beta - 1.0)*c
    for i in range(len(ag)-1):
        dp = -(ag[i+1] - ag[i]) * m
        dpe = dp + a * vel[i] + b*acc[i]

        u_trial = dis[i]
        if i == 1:
            kt = k
        else:
            force_p = force[i-1]
            statev, kt, force_p = Bilinear_Ki(props, force_p, dis[i-1], du, kt, statev)
        ke = kt + gama/beta/dt*c + m/beta/dt/dt
        Res = dpe
        du = 0
        force_c = force[i]
        for j in range(2000):
            du_trial = Res/ke
            du = du + du_trial
            if du_trial / du <= 0.0005:
                break
            force_p = force_c
            statev, kt, force_c = Bilinear_Ki(props, force_c, u_trial, du_trial, kt, statev)
            u_trial = u_trial + du_trial
            ke = kt + gama/beta/dt*c + m/beta/dt/dt
            df = force_c - force_p + (ke-kt) * du_trial
            Res = Res - df
            if j == 1999:
                print(str(T) + '迭代失败')
                return

        dv = gama/beta/dt*du - vel[i]*gama/beta - acc[i]*dt*(1 - gama / 2 / beta)
        da = du/beta/dt/dt - vel[i]/beta/dt - acc[i]/2/beta
        if abs(dis[i] + du) > Dmax:
            Dmax = abs(dis[i] + du)
        if abs(vel[i] + dv) > Vmax:
            Vmax = abs(vel[i] + dv)
        if abs(acc[i] + ag[i+1] + da) > Amax:
            Amax = abs(acc[i] + ag[i+1] + da)
        dis.append(dis[i] + du)
        vel.append(vel[i] + dv)
        acc.append(acc[i] + da)
        force.append(force_c)
        abs_acc.append(acc[i+1] + ag[i+1])
    #print(dis)
    #print(max([abs(x) for x in dis]))
    #print(strength/k)
    #dm = max([abs(x) for x in dis])
    #dy = strength / k
    #print(dm / dy)
    return {'acc': acc,
            'vel': vel,
            'dis': dis,
            'abs_acc': abs_acc,
            'k': k,
            'Dmax': Dmax,
            'Vmax': Vmax,
            'Ab_Amax': Amax,
            'fy': strength,
            'mu': Dmax/(strength/k),
            }


if __name__ =='__main__':
    info = read_record('RSN1_HELENA.A_A-HMC180.AT2')
    res = Nonlinear_SDOF_Analysis(info, 0.5, 0.05, 0.15, 0.01)
    with open('1.txt', 'w') as f:
        for each in res['dis']:
            f.write('%-8.4f\n' % each)
