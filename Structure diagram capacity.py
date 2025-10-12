# Khezr Mohammadaamini
# -----------------------------------------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Input parameters
# -----------------------------------------------------------
A = 4.748e-3       # m²  (Cross-sectional area)
Sigma_max = 500e6  # Pa  (Compressive strength)
Sigma_min = -500e6 # Pa  (Tensile strength)
FS = 1.5           # Factor of Safety
I = 2.231e-5       # m⁴  (Moment of inertia)
t = 0.160          # m   (Steel set height)

# -----------------------------------------------------------
# Derived parameters
# -----------------------------------------------------------

N_max = (A * Sigma_max) / FS
N_min = (A * Sigma_min) / FS
Ncr   = (A * (Sigma_max + Sigma_min)) / (2 * FS)
M_max = ((Sigma_max - Sigma_min) / FS) * (I / t)
M_min = -((Sigma_max - Sigma_min) / FS) * (I / t)

print("N_max =", N_max)
print("N_min =", N_min)
print("M_max =", M_max)
print("M_min =", M_min)
print("Ncr   =", Ncr)

# -----------------------------------------------------------
# Generate Q–N tables (similar to FISH loops)
# -----------------------------------------------------------

# For compressive range (0 → N_max)
N1 = np.linspace(0, N_max, 1000)
Q1 = np.sqrt(np.maximum(0, (4.0 * Sigma_max * A * Sigma_max * A - 4.0 * N1 * FS * Sigma_max * A) / (9.0 * FS**2)))
Q1_neg = -Q1

# For tensile range (N_min → 0)
N2 = np.linspace(N_min, 0, 1000)
Q2 = np.sqrt(np.maximum(0, (4.0 * Sigma_min * A * Sigma_min * A - 4.0 * N2 * FS * Sigma_min * A) / (9.0 * FS**2)))
Q2_neg = -Q2

# Combine and save
Qn_data = pd.DataFrame({
    "N": np.concatenate((N2, N1)),
    "Q": np.concatenate((Q2, Q1))
})
Qn_data.to_csv("Qn_table.csv", index=False)
print("✅ Q–N table saved to 'Qn_table.csv'")

# Plot Q–N curves
# -----------------------------------------------------------

plt.figure(figsize=(8,7))
plt.plot(Q1, N1, 'b', label='Q–N (Compression)')
plt.plot(Q1_neg, N1, 'b')
plt.plot(Q2, N2, 'r', label='Q–N (Tension)')
plt.plot(Q2_neg, N2, 'r')

plt.xlabel('Shear Force Q (N)')
plt.ylabel('Axial Force N (N)')
plt.title('Q–N Interaction Diagram')
plt.legend()
plt.grid(True)

dat1 = pd.read_csv("sm-qn.csv", header=None)

x = dat1[0].astype(float).tolist()
y = dat1[1].astype(float).tolist()
plt.plot(x, y, 'r+', label='sm-n.csv data')
plt.legend()
plt.legend(loc='upper right')
plt.show()

# -----------------------------------------------------------
# Plot M–N envelope
# -----------------------------------------------------------
plt.figure(figsize=(8,7))
plt.plot([0, M_min, 0, M_max, 0], [N_max, 0, N_min, 0, N_max], 'r', label='M–N Envelope (Analytic)')

plt.xlabel('Moment M (N·m)')
plt.ylabel('Axial Force N (N)')
plt.title('Moment–Axial (M–N) Interaction')
plt.legend()
plt.grid(True)

dat2 = pd.read_csv("sm_n.csv", header=None)

x = dat2[0].astype(float).tolist()
y = dat2[1].astype(float).tolist()
plt.plot(x, y, 'b+', label='sm-n.csv data')
plt.legend()
plt.legend(loc='upper right')

plt.show()

print("✅ All calculations and plots completed successfully.")
