import numpy as np
import math
import RectPatch
from RectPatch import DesignPatch, PatchFunction, PatchEHPlanePlot, SurfacePlot, cart2sph1, sph2cart1, SurfacePlot_dB
import ArrayFactor
from ArrayFactor import CalculateRelativePhase
import matplotlib.pyplot as plt

def FieldSumPatch(ElementArray, Freq, W, L, h, Er):
    """
    Summation of field contributions from each patch element in array, at frequency freq for theta 0°-95°, phi 0°-360°.
    Element = xPos, yPos, zPos, ElementAmplitude, ElementPhaseWeight
    Returns arrayFactor[theta, phi, elementSum]
    """

    arrayFactor = np.ones((360, 90))

    Lambda = 3e8 / Freq

    for theta in range(90):
        for phi in range(360):                                                                                                      # For all theta/phi positions
            elementSum = 1e-9 + 0j

            xff, yff, zff = sph2cart1(999, math.radians(theta), math.radians(phi))                                                  # Find point in far field

            for element in ElementArray:                                                                                            # For each element in array, find local theta/phi, calculate field contribution and add to summation for point
                xlocal = xff - element[0]
                ylocal = yff - element[1]                                                                                           # Calculate local position in cartesian
                zlocal = zff - element[2]

                r, thetaLocal, phiLocal = cart2sph1(xlocal, ylocal, zlocal)                                                         # Convert local position to spherical

                patchFunction = PatchFunction(math.degrees(thetaLocal), math.degrees(phiLocal), Freq, W, L, h, Er)            # Patch element pattern for local theta, phi

                if patchFunction != 0:                                                                                              # Sum each elements contribution
                    relativePhase = CalculateRelativePhase(element, Lambda, math.radians(theta), math.radians(phi))                 # Find relative phase for current element
                    elementSum += element[3] * patchFunction * math.e ** ((relativePhase + element[4]) * 1j)                        # Element contribution = Amp * e^j(Phase + Phase Weight)

            arrayFactor[phi][theta] = elementSum.real

    return arrayFactor


if __name__ == "__main__":
    
    freq = 14e9
    W = 10.7e-3
    L = 10.47e-3
    h = 3e-3
    Er = 2.5
    
    #ElementArray = np.array([[0,0,0,1,0], [60e-3,0,0,1,0], [60e-3,60e-3,0,1,0], [0,60e-3,0,1,0]])
    ElementArray = np.genfromtxt("ArrayLayouts/array-10x1.dat")
    fsp = FieldSumPatch(ElementArray, freq, W, L, h, Er)
    #np.savetxt('af.csv', af, delimiter=',')
    #afdb = 20*np.log10(af)

    SurfacePlot_dB(20 * np.log10(abs(fsp)), freq, 0, 0, 0, 0)

    Xtheta = np.linspace(0, 90, 90)
    #plt.plot(Xtheta, af[90, :], label="H-plane (Phi=90°)")
    #plt.plot(Xtheta, af[0, :], label="E-plane (Phi=0°)")
    plt.plot(Xtheta, 20 * np.log10(abs(fsp[90, :])), label="H-plane (Phi=90°)")                          # Log = 20 * log10(E-field)
    plt.plot(Xtheta, 20 * np.log10(abs(fsp[0, :])), label="E-plane (Phi=0°)")
    plt.ylabel('Array Pattern (dB)')
    plt.xlabel('Theta (degs)')                                                                                  # Plot formatting
    #plt.title("Patch: \nW=" + str(W) + " \nL=" + str(L) +  "\nEr=" + str(Er) + " h=" + str(h) + " \n@" + str(Freq) + "Hz")
    #plt.ylim(-40)
    #plt.xlim((0, 90))
    #start, end = plt.xlim()
    #plt.xticks(np.arange(start, end, 5))
    #plt.grid(b=True, which='major')
    plt.legend()
    plt.show()

    #np.savetxt('theta.csv', Xtheta, delimiter=',')
    #np.savetxt('af-python.csv', 20 * np.log10(abs(af[0, :])) -20, delimiter=',')