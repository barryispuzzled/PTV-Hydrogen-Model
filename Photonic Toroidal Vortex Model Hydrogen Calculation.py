import math
from datetime import datetime

def run_ptv_calculation():
    # --- PHYSICAL CONSTANTS (CODATA June 2019) ---
    alpha1 = 0.0072973525693
    alpha2 = alpha1 * alpha1
    Rydfreq = 3289841960.2508
    Me = 9.1093837015 * 10**-31
    clight = 299792458
    M = 0.000544617021488
    pi = 3.141592653589793
    h = 6.62607015 * 10**-34

    # Experimental NIST Data
    states = {
        1: ("1S1/2", 3288087922.4160, 3288086502.0102, 1, 0, 0.5),
        2: ("2S1/2", 822025577.0922, 822025399.5354, 2, 0, 0.5),
        3: ("3S1/2", 365343617.9043, 365343565.2949, 3, 0, 0.5),
        4: ("4S1/2", 205505309.9525, 205505287.7579, 4, 0, 0.5),
        5: ("5S1/2", 131523180.9882, 131523169.6246, 5, 0, 0.5),
        6: ("6S1/2", 91335431.6017, 91335425.0256, 6, 0, 0.5),
        7: ("2P1/2", 822026546.1359, 822026486.9664, 2, 1, -0.5),
        8: ("3P1/2", 365343906.4710, 365343888.9392, 3, 1, -0.5),
        9: ("4P1/2", 205505431.9299, 205505424.5337, 4, 1, -0.5),
        10: ("5P1/2", 131523243.5005, 131523239.7136, 5, 1, -0.5),
        11: ("6P1/2", 91335467.7972, 91335465.6058, 6, 1, -0.5),
        12: ("7P1/2", 67103543.4423, 67103542.0628, 7, 1, -0.5),
        13: ("2P3/2", 822015547.4967, 822015523.8451, 2, 1, 0.5),
        14: ("3P3/2", 365340647.6119, 365340640.6040, 3, 1, 0.5),
        15: ("4P3/2", 205504057.1004, 205504054.1439, 4, 1, 0.5),
        16: ("5P3/2", 131522539.5886, 131522538.0749, 5, 1, 0.5),
        17: ("6P3/2", 91335060.4413, 91335059.5653, 6, 1, 0.5),
        18: ("7P3/2", 67103286.9157, 67103286.3640, 7, 2, -0.5),
        19: ("3D3/2", 365340651.1931, 365340646.9865, 3, 2, -0.5),
        20: ("4D3/2", 205504058.6495, 205504056.8748, 4, 2, -0.5),
        21: ("5D3/2", 131522540.3924, 131522539.4837, 5, 2, -0.5),
        22: ("6D3/2", 91335060.9101, 91335060.3843, 6, 2, -0.5),
        23: ("7D3/2", 67103287.2124, 67103286.8813, 7, 2, -0.5),
        24: ("8D3/2", 51375939.6283, 51375939.4065, 8, 2, -0.5),
        25: ("3D5/2", 365339566.8037, 365339564.1003, 3, 2, 0.5),
        26: ("4D5/2", 205503601.1722, 205503600.0315, 4, 2, 0.5),
        27: ("5D5/2", 131522306.1639, 131522305.5800, 5, 2, 0.5),
        28: ("6D5/2", 91334925.3612, 91334925.0233, 6, 2, 0.5),
        29: ("7D5/2", 67103201.8522, 67103201.6394, 7, 2, 0.5),
        30: ("8D5/2", 51375882.4437, 51375882.3011, 8, 2, 0.5),
    }

    output_file = "PTV.txt"
    print(f"Writing results to {output_file}...")

    with open(output_file, "w") as ptv:
        ptv.write(f"\r\n{datetime.now().strftime('%H:%M:%S  %b %d, %Y')}\r\n\r\n")
        ptv.write("*** A program to calculate hydrogen atom hyperfine levels using the new photonic toroidal vortex (PTV) model ***\r\n")
        ptv.write("              written in Python based on logic by Dr Barry R. Clarke\r\n\r\n")
        ptv.write(" All equations and tables refer to the following papers.\r\n\r\n")
        ptv.write(" PTV1: Clarke, Barry R., A photonic toroidal vortex model of the hydrogen atom fine structure, Quantum Studies:\r\n")
        ptv.write("       Mathematics and Foundations, 12, article 18 (2025)\r\n")
        ptv.write(" PTV2: Clarke, Barry R., The hydrogen atom hyperfine structure from a photonic toroidal vortex model,\r\n")
        ptv.write("       under peer review (2026)\r\n\r\n")
        ptv.write(" EXPERIMENTAL DATA\r\n")
        ptv.write(" Horbatsch, M. and E. A. Hessels, Tabulation of the bound-state energies of atomic hydrogen, Physical Review A 93, 022513 (2016)\r\n\r\n")
        ptv.write("_____________________________________________________________________________________________________\r\n\r\n")

        for nn in range(1, 31):
            label, nist0, nist1, n, L, Sp = states[nn]

            # Hyperfine multipliers
            if 1 <= nn <= 6: mult = 1
            elif 7 <= nn <= 12: mult = 1/3
            elif 13 <= nn <= 18: mult = 2/15
            elif 19 <= nn <= 24: mult = 2/25
            else: mult = 25/486

            na = L + Sp + 0.5
            nr = n - na
            Y = nr + math.sqrt(na**2 - alpha2)
            Y2 = Y * Y

            # Brute force coordinate search logic from BASIC
            x_store = 0.0
            for acc in range(-1, 15):
                store_mark = None
                for msearch in range(0, 151):
                    x_inc = msearch / (10**acc)
                    x_trial = x_store + x_inc
                    r2 = alpha1 / math.sqrt(1 - alpha2)
                    fx_term = n * math.sqrt(1 + alpha2 / Y2) * (1 - alpha2)
                    fint1 = 1 - M
                    fint2 = 1 + (x_trial**2 + r2**2) / (fint1**2)
                    fint3 = 1 / (fint1 * math.sqrt(fint2))
                    fint4 = r2**2 * (x_trial**2 + fint1**2) / (fint1**4 * fint2**2)
                    fint = fint3 * (1 + (0.75 * fint4) + (1.640625 * (fint4**2)))
                    fx = (1 / fx_term) * fint
                    felectron = (0.5 / (Y2 + alpha2)) * (1 / math.sqrt(1 - alpha2 / (2 * (Y2 + alpha2))))
                    mark = 0 if (fx - felectron) < 0 else 1
                    if msearch > 0 and mark != store_mark:
                        x_store = x_store + (msearch - 1) / (10**acc)
                        break
                    store_mark = mark
            
            dist = math.sqrt(x_store**2 + (1 - M)**2)

            # Physics Calculations
            redmass1 = 1 / (1 + M)
            som_energy = (1 - (1 / math.sqrt(1 + (alpha1 / (math.sqrt(na**2 - alpha2) + nr))**2))) * redmass1 * Rydfreq * 2 / alpha2
            evel3 = (1 + M/n)**2 if 1 <= nn <= 6 else (1 + M/(2*n))**2
            e3e = math.sqrt(1 - (alpha2 / Y2) / (2 * (1 + alpha2 / Y2)))
            v_energy0 = Rydfreq / (e3e * (Y2 + alpha2))
            e3ee_arg = (alpha2 * evel3) / (nr + math.sqrt(na**2 - alpha2 * evel3))**2
            e3ee = math.sqrt(1 - e3ee_arg / (2 * (1 + e3ee_arg)))
            redmass3 = 1 / (1 + M / e3ee)
            v_energy1 = (redmass3 * Rydfreq) / (((nr + math.sqrt(na**2 - alpha2 * evel3))**2 + alpha2 * evel3) * e3ee)
            v_energy2 = v_energy1 * evel3
            hf_shift = mult * evel3 * redmass3 * ((10**-6 / h) * Me * alpha1**5 * n**2 * pi * clight**2) / ((dist**3) * (1 - alpha2)**1.5 * math.sqrt(1 + alpha2 / Y2) * Y2 * math.sqrt(2))

            if 1 <= nn <= 6: exitsh1 = 28724616.751299 * (dist**-3.0000402386587108)
            elif 7 <= nn <= 12: exitsh1 = 14326679.482869 * (dist**-3.000001713270)
            elif 13 <= nn <= 18: exitsh1 = 14329441.683161 * (dist**-3.0000113254009108)
            elif 19 <= nn <= 24: exitsh1 = 14328449.553029 * (dist**-3.0000146529976112)
            else: exitsh1 = 14329157.347521 * (dist**-3.0000169617890916)

            # --- OUTPUT FORMATTING MIRRORING PTV.TXT ---
            ptv.write(f"{label}\r\n")
            ptv.write("QUANTUM NUMBERS\r\n")
            ptv.write(f"(1) Quantum Mechanics: n = {int(n)} L = {int(L)} Sp = {Sp}   Sommerfeld: n(phi) = {na} n(r) = {nr}\r\n\r\n")
            
            ptv.write("BOUND STATE DISTANCES (as a muliplier of a fixed radius for all states)\r\n")
            ptv.write(f"(2) Horizontal distance between p-e Sp-2 centers  = {x_store:19.16f}   PTV1 x-bar Eqs (84)/(86)\r\n")
            ptv.write(f"(3) Shortest distance between p-e Sp-2 centers    = {dist:19.16f}   PTV1 d2-bar Eq (86)\r\n\r\n")
            
            ptv.write("FINE STRUCTURE (ALL ENERGIES IN MHz)\r\n")
            ptv.write(f"(4) Sommerfeld f-s energy no adjustments     = {som_energy/redmass1:19.12f}\r\n")
            ptv.write(f"(5) New PTV f-s energy no adjustments        = {v_energy0:19.12f}   PTV1 Eq (54)\r\n")
            ptv.write(f"(6) Sommerfeld f-s energy reduced mass        = {som_energy:19.12f}\r\n")
            ptv.write(f"(7) Relativistic reduced mass multiplier     = {redmass3:19.16f}   PTV2 (1 + Mr)^(-1) Eq (14)\r\n")
            ptv.write(f"(8) New PTV f-s energy relativ. reduced mass = {v_energy1:19.12f}   PTV2 (D = 1) Eq (14)\r\n")
            ptv.write(f"(9) Increased velocity D multiplier          = {evel3:19.16f}   PTV2 Eq (15)\r\n")
            ptv.write("(10) New PTV f-s energy with relativistic reduced mass \r\n")
            ptv.write(f"     and increased velocity                  = {v_energy2:19.12f}   PTV2 Eqs (14)/(15)\r\n\r\n")
            
            ptv.write("EXIT POTENTIAL RED-SHIFT (ALL ENERGIES IN MHz)\r\n")
            ptv.write(f"(11) Average h-f energy                        = {((nist0 + nist1)/2):19.12f}   Horbatsch and Hessels (2016)\r\n")
            ptv.write("(12) Difference between PTV f-s with corrections\r\n")
            ptv.write(f"     and average h-f energy                    = {v_energy2 - (nist0 + nist1)/2:19.12f}\r\n")
            ptv.write(f"(13) Exit potential red-shift                  = {exitsh1:19.12f}   PTV2 Eq (24) & Table 2\r\n")
            ptv.write(f"(14) Adjusted f-s frequency after exit shift   = {v_energy2 - exitsh1:19.12f}   PTV2 Eqs (14)/(15) minus (24)\r\n\r\n")
            
            ptv.write("HYPERFINE SHIFT (ALL ENERGIES IN MHz)\r\n")
            ptv.write(f"(15) Upper h-f energy from experiment      = {nist0:19.12f}       Horbatsch and Hessels (2016)\r\n")
            ptv.write(f"(16) Upper h-f energy calculated from PTV  = {v_energy2 + hf_shift - exitsh1:19.12f}       PTV2 Eqs (14)/(15) minus (24) plus (36)\r\n")
            ptv.write(f"(17) Deviation of PTV from experiment is one part in {abs(nist0/(nist0 - (v_energy2 + hf_shift - exitsh1))):22.12f}\r\n")
            ptv.write(f"(18) Lower h-f energy from experiment      = {nist1:19.12f}       Horbatsch and Hessels (2016)\r\n")
            ptv.write(f"(19) Lower h-f energy calculated from PTV  = {v_energy2 - hf_shift - exitsh1:19.12f}       PTV2 Eqs (14)/(15) minus (24) minus (36)\r\n")
            ptv.write(f"(20) Deviation of PTV from experiment is one part in {abs(nist1/(nist1 - (v_energy2 - hf_shift - exitsh1))):22.12f}\r\n")
            ptv.write(f"(21) H-f difference taken from experiment  = {nist0 - nist1:19.12f}\r\n")
            ptv.write(f"(22) H-f difference calculated from PTV    = {2 * hf_shift:19.12f}       PTV2 twice Eq (36)\r\n\r\n")
            ptv.write("_____________________________________________________________________________________________________\r\n\r\n")

if __name__ == "__main__":
    run_ptv_calculation()
