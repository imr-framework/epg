# Implements 1-look-ahead algorithm for EPG and tests it
from epg import Sequence as EPG_Sequence
from epg import grad_shift, rf_rotation, relaxation
import numpy as np
from math import degrees
from epg import EPG


def generate_one_look_ahead_sequence(echoes, tau, t1t2=(0,0)):
    # Implements 1-look-ahead algorithm from Weigel
    # For convenience, always start with a (a0=90, phi=90) pulse
    # Spacing of rf[0] to rf[1] is 1 tau
    # Spacing from rf[i] to rf[i+1] is 2 tau for all i > 0
    #     and follow with n phi=90 pulses of FA {a1, a2, a3, ...} to get real echoes only
    N = len(echoes)
    echoes = -1j*np.array(echoes)
    phi = 0

    # Set up grad, events, and timing

    epg_grads = (2*N+1)*[1]
    epg_events =  ['rf','grad','relax']
    epg_events.extend(N*['rf','grad','relax','grad','relax'])
    epg_timing = np.array([0,tau,tau])
    epg_timing = np.append(epg_timing, tau + np.array(N*[0,tau,tau,2*tau,2*tau] + np.repeat(2*tau*np.arange(N),5)))

    epg_seq = EPG_Sequence(rf=[(phi,90)], grad=epg_grads, events=epg_events, time=epg_timing, t1t2=t1t2)
    omega = np.matrix([[0], [0], [1]])

    E1 = np.exp(-tau/t1t2[0])if t1t2[0] > 0 else 1
    E2 = np.exp(-tau/t1t2[1])if t1t2[1] > 0 else 1

    # First rf pulse
    omega = rf_rotation(phi, 90)*omega

    for n in range(N):
        # Retrieve previous state for calculating flip angle for next echo
        alpha_next = one_look_ahead_update(I_next=echoes[n], state_prev=omega, E1=E1, E2=E2)
        print(alpha_next)
        # Simulate : {grad, relax, rf, grad, relax}
        omega = grad_shift(1, omega)
        omega = relaxation(tau, t1t2[0], t1t2[1], omega)
        omega = rf_rotation(phi, alpha_next)*omega # rf_rotation outputs an np.matrix so * works as matrix multiplication
        omega = grad_shift(1, omega)
        omega = relaxation(tau, t1t2[0], t1t2[1], omega)

        # Add calculated angle to sequence's rf info
        epg_seq.rf.append((phi, alpha_next))

    epg = EPG(seq=epg_seq)


    print(f"Look-ahead generated for {N} pulses")
    return epg


def one_look_ahead_update(I_next, state_prev, E1, E2):
    #  Retrieve : Zm, Fm, I_prev
    Fm = 0
    Zm = 0
    I_prev = state_prev[0,0]
    if state_prev.shape[1] >= 2:
        Zm = state_prev[2,1]
    if state_prev.shape[1] >= 3:
        Fm = state_prev[1,2]

    a1 = \
        2*np.arctan((-1j*E1*E2*Zm + np.sqrt(-(E1*E2*Zm)**2 - ((E2**2)*Fm - I_next)*((E2**2)*I_prev - I_next)))
                    /((E2**2)*I_prev - I_next))
    a2 = \
        2*np.arctan((-1j*E1*E2*Zm - np.sqrt(-(E1*E2*Zm)**2 - ((E2**2)*Fm - I_next)*((E2**2)*I_prev - I_next)))
                    /((E2**2)*I_prev - I_next))

    print(f"The originals: a1={a1} a2={a2}")

    a1 = degrees(a1)
    a2 = degrees(a2)
    if abs(a1) < 1e-12 : a1 = a2
    if abs(a2) < 1e-12 : a2 = a1

    a = a1 if (abs(a1) > abs(a2)) else a2
    return a


if __name__ == '__main__':
    myepg = generate_one_look_ahead_sequence(echoes=6*[0.2], tau=20, t1t2=(0,0))
    myepg.simulate()
    myepg.display()
    echoes = myepg.plot_echoes()