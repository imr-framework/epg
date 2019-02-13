
# Pre-defined pulse sequences
# Created by Gehua Tong, Feb 13, 2019

import epg
import numpy as np


def make_gre_sequence(alpha,N,TR,t1t2=(0,0)):
    rf = N*[(0,alpha)]
    grad = N*[1]
    events = N*['rf','grad','relax']
    time = N*[0,TR,TR] + np.repeat(TR*np.arange(N),3)
    print(time)
    return epg.Sequence(rf,grad,events,time,t1t2,"Gradient Echo")


def make_tse_sequence(alpha,N,esp,use_y90=1,t1t2=(0,0)):
    events = ['rf','grad','relax']
    events.extend((N-1)*['rf','grad','relax','grad','relax'])
    time = np.array([0,esp/2,esp/2])
    time = np.append(time, esp/2 + np.array((N-1)*[0,esp/2,esp/2,esp,esp] + np.repeat(esp*np.arange(N-1),5)))
    grad = (2*N-1)*[1]
    rf = []
    if use_y90:
        rf.append((90,90))
    else:
        rf.append((0,alpha))

    rf.extend((N-1)*[(0,alpha)])
    return epg.Sequence(rf,grad,events,time,t1t2,"Turbo Spin Echo")


def make_he_sequence(phis,alphas,esp,t1t2=(0,0)):
    if len(phis) != len(alphas):
        raise Exception("phis and alphas must have the same length")
    rf_var_p = []
    rf_var_n = []
    m = len(phis)
    for k in range(m):
        rf_var_p.append((phis[k],alphas[k]))
        rf_var_n.append((-phis[-k-1],-alphas[-k-1]))
    rf = [(90,90)]
    rf.extend(rf_var_p)
    rf.append((0,180))
    rf.extend(rf_var_n)
    grad = (4*m+2)*[1]
    events = ['rf','grad','relax']
    events.extend((2*m)*['rf','grad','relax','grad','relax'])
    events.extend(['rf','grad','relax'])

    time = np.array([0,esp/2,esp/2])
    time = np.append(time, esp/2 + np.array(2*m*[0,esp/2,esp/2,esp,esp]+np.repeat(esp*np.arange(2*m),5)))
    time = np.append(time,[time[-1],time[-1]+esp/2,time[-1]+esp/2])

    return epg.Sequence(rf,grad,events,time,t1t2,"Hyperecho")


def make_star_sequence(alphas,dt,t1t2=(0,0)):
    # Note: gradients are zero
    rf = []
    m = len(alphas)
    for k in range(m):
        rf.append((0,alphas[k]))

    grad = m*[0]
    events = m*['rf','grad','relax']
    time = m*[0,dt,dt]+np.repeat(dt*np.arange(m),3)
    return epg.Sequence(rf,grad,events,time,t1t2,"STAR")
