#Lined Circular Tunnel
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Parameters
def responses(no_slip=True):
    R = 5.0
    K = 0.5
    P = 6e5
    bulk = 5e7
    shear = 1.791e7
    E = 9.0 * bulk * shear / (3.0 * bulk + shear)
    nu = 0.5 * (3.0 * bulk - 2.0 * shear) / (3.0 * bulk + shear)
    Es = 2.5e10
    nus = 0.15
    As = 0.125
    Is = As**3 / 12.0
    Cs = E * R * (1.0 - nus**2) / (Es * As * (1.0 - nu**2))
    Fs = E * R**3 * (1.0 - nus**2) / (Es * Is * (1.0 - nu**2))

    if no_slip:
        Beta = (6.0 + Fs) * Cs * (1.0 - nu) + 2.0 * Fs * nu
        Beta = Beta / (3.0 * Fs + 3.0 * Cs + 2.0 * Cs * Fs * (1.0 - nu))
        bs2 = 2.0 * Cs * (1.0 - nu) + 8.0 * nu - 12.0 * Beta - 6.0 * Beta * Cs * (1.0 - nu)
        bs2 = Cs * (1.0 - nu) / bs2
        as0 = Cs * Fs * (1.0 - nu) / (Cs + Fs + Cs * Fs * (1.0 - nu))
        as2 = Beta * bs2
    else:
        as0 = Cs * Fs * (1.0 - nu) / (Cs + Fs + Cs * Fs * (1.0 - nu))
        as2 = (Fs + 6.0) * (1.0 - nu) / (2.0 * Fs * (1.0 - nu) + 6.0 * (5.0 - 6.0 * nu))

    if no_slip:
        U1 = 0.5 * (1.0 + K) * as0
        U2 = 0.5 * (1.0 - K) * (4.0 * (1.0 - nu) * bs2 - 2.0 * as2)
        V2 = (-1.0 + K) * (as2 + (1.0 - 2.0 * nu) * bs2)
        T1 = (1.0 + K) * (1.0 - as0)
        T2 = (1.0 - K) * (1.0 + 2.0 * as2)
        M2 = 0.5 * (1.0 - K) * (1.0 - 2.0 * as2 + 2.0 * bs2)
        sigR1 = (1.0 + K) * (1.0 - as0)
        sigR2 = (-1.0 + K) * (1.0 - 6.0 * as2 + 4.0 * bs2)
        tauRT2 = (1.0 - K) * (1.0 + 6.0 * as2 - 2.0 * bs2)
    else:
        U1 = 0.5 * (1.0 + K) * as0
        U2 = (-1.0 + K) * ((5.0 - 6.0 * nu) * as2 - (1.0 - nu))
        V2 = 0.5 * (1.0 - K) * ((5.0 - 6.0 * nu) * as2 - (1.0 - nu))
        T1 = (1.0 + K) * (1.0 - as0)
        T2 = (1.0 + K) * (1.0 - 2.0 * as2)
        M2 = (1.0 - K) * (1.0 - 2.0 * as2)
        sigR1 = (1.0 + K) * (1.0 - as0)
        sigR2 = (-1.0 + K) * (3.0 - 6.0 * as2)
        tauRT2 = 0.0

    theta_deg = np.linspace(0, 90, 200)
    theta = np.deg2rad(theta_deg)

    ua = (P * R * (1.0 + nu) / E) * (U1 + U2 * np.cos(2.0 * theta))
    va = (P * R * (1.0 + nu) / E) * (V2 * np.sin(2.0 * theta))

    sigRa = 0.5 * P * (sigR1 + sigR2 * np.cos(2.0 * theta))
    tauRTa = 0.5 * P * tauRT2 * np.sin(2.0 * theta)

    Ta = 0.5 * P * R * (T1 + T2 * np.cos(2.0 * theta))
    Ma = 0.5 * P * R**2 * M2 * np.cos(2.0 * theta)

    df = pd.DataFrame({
        'theta_deg': theta_deg,
        'ua': ua,
        'va': va,
        'sigR': sigRa,
        'tauRT': tauRTa,
        'T': Ta,
        'M': Ma
    })

    df.to_csv('analytic_results_0_90.csv', index=False)

    # --- Plot results ---
    plt.figure(figsize=(12, 8))
    plt.subplot(2, 2, 1)
    plt.plot(theta_deg, ua, label='ua')
    plt.plot(theta_deg, va, label='va')
    plt.xlabel('Theta (deg)')
    plt.ylabel('Displacement (m)')
    plt.title('Support Displacements')
    plt.legend()

    plt.subplot(2, 2, 2)
    plt.plot(theta_deg, sigRa, label='sigR')
    plt.plot(theta_deg, tauRTa, label='tauRT')
    plt.xlabel('Theta (deg)')
    plt.ylabel('Stress (Pa)')
    plt.title('Interface Stresses')
    plt.legend()

    plt.subplot(2, 2, 3)
    plt.plot(theta_deg, Ta, label='T')
    plt.xlabel('Theta (deg)')
    plt.ylabel('Axial Force (N/m)')
    plt.title('Stress Resultant T')
    plt.legend()

    plt.subplot(2, 2, 4)
    plt.plot(theta_deg, Ma, label='M')
    plt.xlabel('Theta (deg)')
    plt.ylabel('Moment (N·m/m)')
    plt.title('Stress Resultant M')
    plt.legend()

    plt.tight_layout()
    plt.show()

    print('Results saved to analytic_results_0_90.csv')

responses(no_slip=True)
