import numpy as np
import math
import RectPatch
from RectPatch import SurfacePlot_dB
import matplotlib.pyplot as plt

def ArrayFactor(ElementArray, Freq):
    """
    Summation of field contributions from each element in array, at frequency freq at theta 0°-95°, phi 0°-360°.
    Element = xPos, yPos, zPos, ElementAmplitude, ElementPhaseWeight
    Returns arrayFactor[theta, phi, elementSum]
    """

    arrayFactor = np.ones((360, 90))

    Lambda = 3e8 / Freq

    for theta in range(90):
        for phi in range(360):                                                                                                      # For all theta/phi positions
            elementSum = 1e-9 + 0j

            for element in ElementArray:                                                                                            # Summation of each elements contribution at theta/phi position.
                relativePhase = CalculateRelativePhase(element, Lambda, math.radians(theta), math.radians(phi))                     # Find relative phase for current element
                elementSum += element[3] * math.e ** ((relativePhase + element[4]) * 1j)                                            # Element contribution = Amp * e^j(Phase + Phase Weight)

            arrayFactor[phi][theta] = abs(elementSum)
            #arrayFactor[phi][theta] = elementSum.real

    return arrayFactor


def CalculateRelativePhase(Element, Lambda, theta, phi):
    """
    Incident wave treated as plane wave. Phase at element is referred to phase of plane wave at origin.
    Element = xPos, yPos, zPos, ElementAmplitude, ElementPhaseWeight
    theta & phi in radians
    See Eqn 3.1 @ https://theses.lib.vt.edu/theses/available/etd-04262000-15330030/unrestricted/ch3.pdf
    """
    phaseConstant = (2 * math.pi / Lambda)

    xVector = Element[0] * math.sin(theta) * math.cos(phi)
    yVector = Element[1] * math.sin(theta) * math.sin(phi)
    zVector = Element[2] * math.cos(theta)

    phaseOfIncidentWaveAtElement = phaseConstant * (xVector + yVector + zVector)

    return phaseOfIncidentWaveAtElement


if __name__ == "__main__":
    #ElementArray = np.array([[0,0,0,1,0], [60e-3,0,0,1,0], [60e-3,60e-3,0,1,0], [0,60e-3,0,1,0]])
    ElementArray = np.genfromtxt("ArrayLayouts/array-10x1.dat")
    af = ArrayFactor(ElementArray, 2.4e9)
    #np.savetxt('af.csv', af, delimiter=',')
    #afdb = 20*np.log10(af)

    SurfacePlot_dB(20 * np.log10(abs(af)), 2.4e9, 0, 0, 0, 0)

    Xtheta = np.linspace(0, 90, 90)
    #plt.plot(Xtheta, af[90, :], label="H-plane (Phi=90°)")
    #plt.plot(Xtheta, af[0, :], label="E-plane (Phi=0°)")
    plt.plot(Xtheta, 20 * np.log10(abs(af[90, :])), label="H-plane (Phi=90°)")                          # Log = 20 * log10(E-field)
    plt.plot(Xtheta, 20 * np.log10(abs(af[0, :])), label="E-plane (Phi=0°)")
    plt.ylabel('AF (dB)')
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
