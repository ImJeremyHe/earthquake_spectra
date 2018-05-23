"""
参考曲哲老师的Fortran程序：力—变形关系为双折线模型
props—SDOF体系的特性
s—
e—
de—
Et—
"""


def Bilinear_Ki(props, s, e, de, Et, statev):
    E0 = props[0]   #初始刚度
    sy = props[1]   #屈服力
    eta = props[2]  #Strain hardening ratio

    s = s + E0 * de
    Et = E0
    if de >= 0:
        evs = sy + (e + de - sy/E0) * eta * E0
        evE = eta * E0
        if s >= evs:
            s = evs
            Et = evE
    elif de < 0:
        evs = -sy + (e + de + sy/E0) * eta * E0
        evE = eta * E0
        if s <= evs:
            s = evs
            Et = evE
    return statev, Et, s



'''Clough退化双线型模型'''
def Bilinear_Degegrading(props, s, e, de, Et, statev):
    E0 = props[0]
    sy = props[1]
    eta = props[2]

    [emax, emin, ert, srt, erc, src, kon] = statev
    kon = round(kon)

    if kon ==0:
        emax = sy/E0
        emin = -sy/E0
        if de<= 0:
            kon = 1
        else:
            kon = 2
    elif kon == 1 and de < 0:
        kon = 2
        if s > 0:
            erc = e
            src = s
        if e > emax:
            emax = e
    elif kon == 2 and de > 0:
        kon  = 1
        if s < 0:
            ert = e
            srt = s
        if e < emin:
            emin = e

    s = s + E0 * de
    Et = E0
    if de >= 0:
        evs = sy + (e + de - sy/E0) * eta * E0
        evE = eta * E0
        if s >= evs:
            s = evs
            Et = evE
        smax = max(sy, sy + (emax - sy/E0) * eta * E0)
        sres = 0
        eres = ert - (srt - sres) / E0
        if eres <= emax -smax / E0:
            srel = (e + de - eres) / (emax - eres) * (smax - sres) + sres
            if s > srel:
                s = srel
                Et = (smax - sres) / (emax - eres)
    elif de < 0:
        evs = -sy + (e + de + sy/E0) * eta * E0
        evE = eta * E0
        if s <= evs:
            s = evs
            Et = evE

        smin = min(-sy, -sy + (emin + sy/E0) * eta * E0)
        sres = 0
        eres = erc - (src - sres) / E0
        if eres >= emin - smin / E0 :
            srel = (e + de - eres)/(emin - eres) * (smin - sres) + sres
            if s <= srel :
                s = srel
                Et = (smin - sres) / (emin - eres)
    statev = [emax, emin, ert, srt, erc, src, kon]
    return statev, Et, s


