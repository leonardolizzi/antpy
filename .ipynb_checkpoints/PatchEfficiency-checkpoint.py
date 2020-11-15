"""
Returns the efficiency of a rectangular microstrip patch as a percentage. Based on ArrayCalc calc_patchr_eff.m.
References:
Microstrip Antennas, I.J Bahl and P.Bhartia, Published Atrech House, Page 60
Advances in Microstrip and Printed Antennas", Lee and Chen (Ch5)

Some useful numbers :

                CONDUCTORS                                      DIELECTRICS

        Material             Conductivity S/m           Material         Er     Tand

        Perfect              9.90E+99 (lossless)        FR4_Epoxy        4.4    0.02
        Silver               6.29E+07                   Arlon 25FR       3.43   0.0035
        Copper               5.80E+07                   Arlon AD300      3.00   0.003
        Pure Alumin.         3.77E+07                   Arlon AR1000    10.00   0.0035
        Al. 6063-T832        3.08E+07                   Rogers RO3003    3.00   0.0013
        Al. 6061-T6          2.49E+07                   Rogers RO3006    6.15   0.0025
        Brass                1.56E+07                   Rogers RO3010   10.20   0.0035
        Phospor bronze       9.09E+06                   Rogers RO4350    3.48   0.004
        Stainless Steel 302  1.39E+06                   Glass             5.5   0.000
                                                        Plexiglass        3.4   0.001
                                                        Polyamide         4.3   0.004
                                                        Polyester         3.2   0.003
                                                        Polyethylene      2.25  0.001
"""
from math import sqrt, pi
import RectPatch
from RectPatch import DesignPatch, PatchFunction
import Directivity
from Directivity import CalcDirectivity


def CalculatePatchEff(Er, W, L, h, tand, sigma, Freq, VSWR):
    """
    Er: Relative dielectric constant
    W: Patch width (m)
    L: Patch length (m)
    h: dielectric thickness (m)
    tand: Loss tangent of dielectric (units)
    sigma: Conductivity of patch (Siemens/m)
    Freq: Frequency (Hz)
    VSWR: VSWR for bandwidth estimate (ratio). http://www.antenna-theory.com/definitions/vswr.php describes how well the antenna is impedance matched to the line it is connected to.
    """

    if Er <= 1.000001:
        Er = 1.000001

    if tand <= 0.000001:
        tand = 0.000001

    Eo = 8.854185e-12                                                                               # Free space dielectric constant
    Ee = Eo * Er                                                                                    # Effective dielectric constant

    lamba = 3e8 / Freq

    """
    % Calculation for space and surface wave efficiency factor, gives roughly the same results.
    % Reference :  "Advances in Microstrip and Printed Antennas" Lee and Chen(Ch 5)
    % Efficiency due to surface wave component, dominant for larger h/lambda values
    """
    Mur = 1
    n1 = sqrt(Er * Mur)
    ko = 2 * pi * Freq * (sqrt((8.854e-12) * (pi * 4e-7)))
    Lo = lamba

    Psur = (1 / Lo ** 2) * ((ko * h) ** 3) * (60 * pi ** 3 * Mur ** 3 * (1 - 1 / n1 ** 2) ** 3)    # Power radiated as surface wave

    c1 = 1 - 1 / n1 ** 2 + (2 / 5) / n1 ** 4
    Pr = (1 / Lo ** 2) * (ko * h) ** 2 * (80 * pi ** 2 * Mur ** 2 * c1)                                   # Total power radiated

    Effsw = Pr / (Pr + Psur)                                                                        # Efficiency factor for surface wave losses

    """
    % Efficiency due to ohmic and dielectric losses, dominant for smaller h/lambda values
    % ***********************************************************************************
    % Reference : "Microstrip Antennas" Bahl and Bartia
    """
    if W < lamba:
        Rr = 90 * lamba ** 2 / W ** 2                                                               # Radiation resistance for W<lambda
    else:
        Rr = 120 * lamba / W                                                                        # Radiation resistance for W>lambda

    Qr = (3e8 * sqrt(Ee)) / (2 * (Freq / 1e6) * h)                                                    # Quality factor, modified by me, not sure freq was in Ghz, more like MHz !?
    Rc = (1 / sigma) * 0.5 * sqrt(Freq) * (L / W) * Qr ** 2                                         # Equivalent resistance for conductor loss (ohms)
    Rd = (30 * tand / Er) * ((h * lamba) / (L * W)) * Qr ** 2                                       # Equivalent resistance for dielectric loss (ohms)

    Rtot = Rr + Rd + Rc                                                                             # Total resistance (ohms)
    Effcd = Rr / Rtot                                                                               # Efficiency factor for combined dielectric and ohmic losses

    Eff1 = Effsw * Effcd
    Eff = Eff1 * 100                                                                                # Total efficiency including ohmic, dielectric and surface wave losses (percent)

    Qt = Qr * Eff1 / (pi)                                                                           # Ref Balanis p762  ( Qtotal = Qradiated*Efficiency ) Not the pi factor, I added that, seems necassary get sensible results using Qr from above !?

    BW = (VSWR - 1) / (Qt * sqrt(VSWR))

    BWf = BW * Freq / 1e6                                                                           # Bandwidth as a frequency span in MHz

    print("Rectangular patch overall efficency " + str(Eff) + "%")
    print("Surface wave efficiency factor " + str(Effsw))
    print("Ohmic and dielectric efficiency factor " + str(Effcd))
    print("BW=" + str(BWf) + "MHz for VSWR=" + str(VSWR) + " at Fo=" + str(Freq / 1e6) + " MHz")

    return Eff


if __name__ == "__main__":
    """
    Calculates efficiencies for two patches @14GHz, one with FR4 and one RO4350.
    """

    freq = 14e9                     # Same parameters for both patches
    h = 1.524e-3
    VSWR = 2.0
    sigma = 5.8e7                   # Copper

    # FR4_Epoxy
    Er = 4.4
    tand = 0.02

    print("\n\nCalculating for FR4 patch.")
    W, L, h, Er = DesignPatch(Er, h, freq)
    eff = CalculatePatchEff(Er, W, L, h, tand, sigma, freq, VSWR)
    CalcDirectivity(eff, PatchFunction, freq, W, L, h, Er)

    # Rogers RO4350
    print("\n\nCalculating for RO4350 patch.")
    Er = 3.48
    tand = 0.004
    W, L, h, Er = DesignPatch(Er, h, freq)
    eff = CalculatePatchEff(Er, W, L, h, tand, sigma, freq, VSWR)
    CalcDirectivity(eff, PatchFunction, freq, W, L, h, Er)
