
# Pre-defined pulse sequences
# Created by Gehua Tong, Feb 13, 2019

import epg
import numpy as np
import math

def make_gre_sequence(alpha,N,TR,t1t2=(0,0),dph=0):
    """Constructs GRE sequence (repeated RF excitation with constant interval TR)

    Inputs
    ------
    alpha : float
        Nutation angle [degrees]
    N : int
        Total number of RF pulses
    TR : float
        Repetition time [ms]
    t1t2 : tuple, optional
        2-tuple (t1, t2) of relaxation times [ms] to be used
        Default is (0,0), which means no relaxation here
    dph : float, optional
        Phase step for RF spoiling option;
        Default is zero - no RF spoiling

    Returns
    -------
    seq : epg.Sequence
        Constructed GRE sequence

    """
    if dph == 0:
        rf = N*[(0,alpha)]
    else:
        phases=  dph*(1 + np.arange(0,N)*np.arange(1,N+1)/2)
        rf = np.zeros((N,2))
        rf[:,0] = phases
        rf[:,1] = alpha*np.ones(N)


    grad = N*[1]
    events = N*['rf','grad','relax']
    time = N*[0,TR,TR] + np.repeat(TR*np.arange(N),3)
    seq = epg.Sequence(rf,grad,events,time,t1t2,"Gradient Echo")

    return seq



def make_tse_sequence(alpha,N,esp,use_y90=True, t1t2=(0,0)):
    """Constructs TSE sequence (repeated RF excitation with constant interval TR)

    Inputs
    ------
    alpha : float
        Nutation angle [degrees] for repeated pulses (2nd to last)
    N : int
        Total number of RF pulses
    esp : float
        Echo spacing [ms]; first interval is esp/2 and second to last intervals are all esp
    use_y90 : boolean, optional
        Whether to use a 90-deg angle with 90-deg phase (y direction) at the beginning.
        Default is True.
        If False, the first pulse is replaced with an alpha-deg pulse with a 0-deg phase (x direction)
    t1t2 : tuple, optional
        2-tuple (t1, t2) of relaxation times [ms] to be used
        Default is (0,0), which means no relaxation here


    Returns
    -------
    seq : epg.Sequence
        Constructed TSE sequence

    """
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
    """Constructs HE (hyperecho) sequence
    # See Hennig, J., & Scheffler, K. (2001). Hyperechoes. Magnetic Resonance in Medicine, 46(1), 6-12.

    Inputs
    ------
    phis : array_like
        RF phases [degrees] to use in the first half
    alpha : array_like
        Nutation phase [degrees] to use in the first half
    esp : float
        Echo spacing [ms]; spacing between all RF pulses except for the first interval
    t1t2 : tuple, optional
        2-tuple (t1, t2) of relaxation times [ms] to be used
        Default is (0,0), which means no relaxation here

    Returns
    -------
    seq : epg.Sequence
        Constructed HE sequence

    """

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


def make_star_sequence(alpha,N,dt,t1t2=(0,0)):
    """Constructs a discrete approximation of a STAR sequence with constant RF excitation
    # Note: gradients are set to zero, so it is effectively equivalent to time-domain simulation with a single isochromat

    Inputs
    ------
    alpha : float
        Nutation angle [degrees]
    N : int
        Total number of RF pulses
    dt : float
        Spacing between RF pulses [ms]
    t1t2 : tuple, optional
        2-tuple (t1, t2) of relaxation times [ms] to be used
        Default is (0,0), which means no relaxation here


    Returns
    -------
    seq : epg.Sequence
        Constructed TSE sequence

    """


    rf = []
    for k in range(N):
        rf.append((0,alpha))

    grad = N*[0]
    events = N*['rf','grad','relax']
    time = N*[0,dt,dt]+np.repeat(dt*np.arange(N),3)
    return epg.Sequence(rf,grad,events,time,t1t2,"STAR")


def make_star_vfa_sequence(alphas,dt,t1t2=(0,0)):
    """Constructs a discrete approximation of a STAR sequence with B1-varying RF excitation (phase is fixed at 0)
    # Note: gradients are set to zero, so it is effectively equivalent to time-domain simulation with a single isochromat

    Inputs
    ------
    alphas : array_like
        Nutation angles [degrees]; its length is the number of RF pulses
    dt : float
        Spacing between RF pulses [ms]
    t1t2 : tuple, optional
        2-tuple (t1, t2) of relaxation times [ms] to be used
        Default is (0,0), which means no relaxation here

    Returns
    -------
    seq : epg.Sequence
        Constructed TSE sequence

    """
    rf = []
    m = len(alphas)
    for k in range(m):
        rf.append((0,alphas[k]))

    grad = m*[0]
    events = m*['rf','grad','relax']
    time = m*[0,dt,dt]+np.repeat(dt*np.arange(m),3)
    return epg.Sequence(rf,grad,events,time,t1t2,"STAR VFA")


def make_bssfp_sequence(alphas, tr, beta=0, t1t2=(0,0),dk=1):
    """Constructs a bssfp sequence with constant phase shift between neighboring pulses

    Inputs
    ------
    alphas : array_like
        Nutation angles [degrees]; its length is the number of RF pulses
    tr : float
        Repetition time [ms]
    beta : float, optional
        Phase shift [degrees] between pulses. Default is 0.
    t1t2 : tuple, optional
        2-tuple (t1, t2) of relaxation times [ms] to be used
        Default is (0,0), which means no relaxation here
    dk : int, optional
        Configuration state shift [a.u.]
        Default is 1.

    Returns
    -------
    seq : epg.Sequence
        Constructed TSE sequence

    """
    n = 0
    m = len(alphas)
    rf = []
    for k in range(len(alphas)):
        rf.append((n*beta, alphas[k]))
        n += 1
    grad = m*[dk]
    events = m*['rf','grad','relax']
    time = m*[0,tr,tr] + np.repeat(tr*np.arange(m),3)
    return epg.Sequence(rf,grad,events,time,t1t2,"bSSFP")

def make_bssfp_sequence_alt(alpha, N, tr, t1t2=(0,0), dk=1):
    """Make a constant-flip angle bssfp sequence with intermediate sampling (at midpoint of TR)
       Note: the first pulse has flip angle of (-alpha/2) and is followed by: alpha, -alpha, alpha, -alpha ...
             (i.e. there is an implicit phase shift of 180 between pulses)
             Also, the spacing between the 1st and 2nd pulse is (TR/2) while all other spacings are TR
             This arrangement allows a relatively high signal and faster convergence to steady state.

    Inputs
    ------
    alpha : float
        Nutation angle [degrees]
    N : int
        Number of RF pulses
    tr : float
        Repetition time [ms]
    t1t2 : tuple, optional
        2-tuple (t1, t2) of relaxation times [ms] to be used
        Default is (0,0), which means no relaxation here
    dk : int, optional
        Configuration state shift [a.u.]
        Default is 1.

    Returns
    -------
    seq : epg.Sequence
        Constructed TSE sequence


    """
    esp = tr
    events = ['rf','grad','relax']
    events.extend((N-1)*['rf','grad','relax','grad','relax'])
    time = np.array([0,esp/2,esp/2])
    time = np.append(time, esp/2 + np.array((N-1)*[0,esp/2,esp/2,esp,esp] + np.repeat(esp*np.arange(N-1),5)))
    grad = (2*N-1)*[dk]


    rf = [(0,-alpha/2)]
    rf.extend(int((N-1)/2)*[(0,alpha),(0,-alpha)])

    return epg.Sequence(rf, grad, events, time, t1t2, 'bSSFP_alt')

