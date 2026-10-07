
# Με βάση τους όγκους C3H7OH και Η2Ο που έχουν προστεθεί σε κάθε διάλυμα, υπολογίζουμε τα moles κάθε
# συστατικού σε κάθε διάλυμα
# και στη συνέχεια το γραμμομοριακό κλάσμα x1 με βάση την C3H7OH

# Τέλος, σε κοινό διάγραμμα Τxy καταγράφουμε το σημείο που αντιστοιχεί σε έκαστο διάλυμα, ήτοι το σημείο που
#ορίζεται από την θερμοκρασία ζέσεως κάθε διαλύματος και το γραμμομοριακό κλάσμα C3H7OH εκάστου διαλύματος.
#Προφανώς με αυτό τον τρόπο σχεδιάζουμε μόνο την καμπύλη υγρού
#Συνδέουμε μεταξύ τους τα σημεία και εκτιμούμε, είτε μαθηματικά είτε γραφικά, τη θερμοκρασία και τη σύσταση του
#αζεοτρόπου

# Δεδομένα εργαστηρίου:


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import sys

sys.stdout.reconfigure(encoding="utf-8")  # type: ignore # για να τυπώνονται τα ελληνικά στο terminal των Windows

# Φυσικές σταθερές (καθαρά υγρά, ~20 °C)
rho_IPA = 0.786      # g/mL, πυκνότητα ισοπροπανόλης
rho_H2O = 0.998      # g/mL, πυκνότητα νερού
M_IPA = 60.10        # g/mol, μοριακή μάζα C3H7OH
M_H2O = 18.015       # g/mol, μοριακή μάζα H2O

# Διάλυμα 1: αρχικά 100 mL ισοπροπανόλης, προσθήκες νερού
V0_IPA_sol1 = 100.0                              # mL
Isopropanol_added_water = [0, 5, 5, 10, 20, 20]  # mL νερού σε κάθε βήμα
Isopropanol_temps = [80.5, 79.5, 79.5, 79.5, 80, 80.5]  # °C

# Διάλυμα 2: αρχικά 100 mL νερού, προσθήκες ισοπροπανόλης
V0_H2O_sol2 = 100.0                              # mL
Water_added_IPA = [0, 20, 20, 10, 10, 10]        # mL ισοπροπανόλης σε κάθε βήμα
Water_temps = [98, 84.5, 81, 81, 81, 81]         # °C


def mol_from_volume(V, rho, M):
    return V * rho / M


# Αθροιστικοί όγκοι σε κάθε βήμα
V_IPA_sol1 = np.full(len(Isopropanol_added_water), V0_IPA_sol1)
V_H2O_sol1 = np.cumsum(Isopropanol_added_water)

V_IPA_sol2 = np.cumsum(Water_added_IPA)
V_H2O_sol2 = np.full(len(Water_added_IPA), V0_H2O_sol2)

# Moles κάθε συστατικού
n_IPA_sol1 = mol_from_volume(V_IPA_sol1, rho_IPA, M_IPA)
n_H2O_sol1 = mol_from_volume(V_H2O_sol1, rho_H2O, M_H2O)
n_IPA_sol2 = mol_from_volume(V_IPA_sol2, rho_IPA, M_IPA)
n_H2O_sol2 = mol_from_volume(V_H2O_sol2, rho_H2O, M_H2O)

# Γραμμομοριακό κλάσμα ισοπροπανόλης x1
x1_sol1 = n_IPA_sol1 / (n_IPA_sol1 + n_H2O_sol1)
x1_sol2 = n_IPA_sol2 / (n_IPA_sol2 + n_H2O_sol2)

df_sol1 = pd.DataFrame({
    "V_IPA (mL)": V_IPA_sol1,
    "V_H2O (mL)": V_H2O_sol1,
    "n_IPA (mol)": n_IPA_sol1,
    "n_H2O (mol)": n_H2O_sol1,
    "x1_IPA": x1_sol1,
    "T (°C)": Isopropanol_temps,
})
df_sol2 = pd.DataFrame({
    "V_IPA (mL)": V_IPA_sol2,
    "V_H2O (mL)": V_H2O_sol2,
    "n_IPA (mol)": n_IPA_sol2,
    "n_H2O (mol)": n_H2O_sol2,
    "x1_IPA": x1_sol2,
    "T (°C)": Water_temps,
})

pd.set_option("display.float_format", "{:.4f}".format)
print("Διάλυμα 1 (αρχικά ισοπροπανόλη, προσθήκη νερού):")
print(df_sol1.to_string(index=False))
print("\nΔιάλυμα 2 (αρχικά νερό, προσθήκη ισοπροπανόλης):")
print(df_sol2.to_string(index=False))

# Κοινή καμπύλη υγρού: ενώνουμε τα σημεία και των δύο διαλυμάτων ταξινομημένα κατά x1
x_all = np.concatenate([x1_sol1, x1_sol2])
T_all = np.concatenate([Isopropanol_temps, Water_temps])
order = np.argsort(x_all)
x_all, T_all = x_all[order], T_all[order]

# Μαθηματική εκτίμηση αζεοτρόπου: προσαρμογή παραβολής στα σημεία κοντά στο ελάχιστο (T <= 82 °C)
mask = T_all <= 82
a, b, c = np.polyfit(x_all[mask], T_all[mask], 2)
x_az = -b / (2 * a)
T_az = c - b**2 / (4 * a)

# Γραφική εκτίμηση: ελάχιστη μετρημένη θερμοκρασία και εύρος x1 όπου εμφανίζεται
T_min_meas = T_all.min()
x_at_min = x_all[T_all == T_min_meas]

# Μετατροπή σύστασης αζεοτρόπου σε κλάσμα μάζας και % v/v (ιδανική ανάμιξη όγκων)
m_IPA_az = x_az * M_IPA
m_H2O_az = (1 - x_az) * M_H2O
w_az = m_IPA_az / (m_IPA_az + m_H2O_az)
V_IPA_az = m_IPA_az / rho_IPA
V_H2O_az = m_H2O_az / rho_H2O
phi_az = V_IPA_az / (V_IPA_az + V_H2O_az)

# Βιβλιογραφικές τιμές (P = 1 atm)
x_az_lit = 0.68
T_az_lit = 80.3

print("\n" + "=" * 60)
print("ΑΠΟΤΕΛΕΣΜΑΤΑ ΑΖΕΟΤΡΟΠΟΥ ΙΣΟΠΡΟΠΑΝΟΛΗΣ (1) – ΝΕΡΟΥ (2)")
print("=" * 60)
print("Γραφική εκτίμηση (από τα πειραματικά σημεία):")
print(f"  T_min = {T_min_meas:.1f} °C  για x1 = {x_at_min.min():.3f} – {x_at_min.max():.3f}")
print("\nΜαθηματική εκτίμηση (παραβολή στα σημεία με T <= 82 °C):")
print(f"  T = {a:.3f}·x1² + ({b:.3f})·x1 + {c:.3f}")
print(f"  Θερμοκρασία αζεοτρόπου  T_az   = {T_az:.2f} °C")
print(f"  Γραμμομοριακό κλάσμα    x1_az  = {x_az:.3f}   (x2 = {1 - x_az:.3f})")
print(f"  Κλάσμα μάζας IPA        w1_az  = {w_az:.3f}   ({100 * w_az:.1f} % w/w)")
print(f"  Κλάσμα όγκου IPA        φ1_az  = {phi_az:.3f}   ({100 * phi_az:.1f} % v/v)")
print("\nΣύγκριση με βιβλιογραφία (1 atm):")
print(f"  x1_lit = {x_az_lit:.2f}, T_lit = {T_az_lit:.1f} °C")
print(f"  Απόκλιση x1: {100 * (x_az - x_az_lit) / x_az_lit:+.1f} %")
print(f"  Απόκλιση T : {T_az - T_az_lit:+.2f} °C ({100 * (T_az - T_az_lit) / T_az_lit:+.1f} %)")
print("\nΤο αζεότροπο είναι ελαχίστου σημείου ζέσεως (θετική απόκλιση από τον νόμο του Raoult).")
print("=" * 60)

fig, ax = plt.subplots(figsize=(8, 5.5))
ax.plot(x_all, T_all, "-", color="0.6", lw=1.5, zorder=1, label="Καμπύλη υγρού (σύνδεση σημείων)")
ax.plot(x1_sol1, Isopropanol_temps, "o", ms=8, color="#2a6fdb", mec="white", mew=1.5, zorder=3,
        label="Διάλυμα 1 (IPA + προσθήκη H$_2$O)")
ax.plot(x1_sol2, Water_temps, "s", ms=8, color="#e07b00", mec="white", mew=1.5, zorder=3,
        label="Διάλυμα 2 (H$_2$O + προσθήκη IPA)")

x_fit = np.linspace(x_all[mask].min(), x_all[mask].max(), 200)
ax.plot(x_fit, np.polyval([a, b, c], x_fit), "--", color="k", lw=1.2, zorder=2,
        label="Παραβολική προσαρμογή (T ≤ 82 °C)")
ax.plot(x_az, T_az, "*", ms=16, color="#c0392b", mec="white", mew=1, zorder=4,
        label=f"Αζεότροπο: x$_1$ ≈ {x_az:.3f}, T ≈ {T_az:.1f} °C")
ax.axvline(x_az, color="#c0392b", lw=0.8, ls=":", zorder=0)
ax.axhline(T_az, color="#c0392b", lw=0.8, ls=":", zorder=0)

ax.set_xlabel("x$_1$ (γραμμομοριακό κλάσμα C$_3$H$_7$OH)")
ax.set_ylabel("T (°C)")
ax.set_title("Διάγραμμα T–x: σύστημα ισοπροπανόλη (1) – νερό (2)")
ax.set_xlim(0, 1)
ax.grid(True, color="0.9")
ax.legend(frameon=False)
fig.tight_layout()
plt.show()



