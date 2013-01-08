from pylab import *
import gkcStyle

from iCode import *

Nx = 513
stdSetup_swITG = { 'Ln' : 1. , 'Ls' : 1./0.2, 'Lx' : 20., 'Debye2' : 0., 'species' : [ Species(m=0.,n=1.,q=-1,name="Adiab"),  Species(n=1.,m=1.,q=1.,eta=5.,name="Ion") ] }

disp = solveDispersion(Nx)
w0 = -0.25 + 0.12j
"""
# std-ITG (1)
w_0, kx, phi_k, x, phi_0  = disp.solveEquation((w0, w0+0.0005, w0+0.0002j), 0.49, stdSetup_swITG)

# sw-ITG (1)
w_1, kx, phi_k, x, phi_1  = disp.solveEquation((w0, w0+0.0005, w0+0.0002j), 2.5, stdSetup_swITG)

# sw-ITG (2)
w_2, kx, phi_k, x, phi_2  = disp.solveEquation((w0, w0+0.0005, w0+0.0002j), 9., stdSetup_swITG)
"""

# sw-ITG (2)
w_2, kx, phi_k, x, phi_2  = disp.solveEquation((w0, w0+0.0005, w0+0.0002j), 0.5, stdSetup_swITG)


########### Plotting #######333
gkcStyle.newFigure(ratio="2.33:1", basesize=10.)

subplot(121)
plot(x, real(phi_0), label="%.3f %.3f" % (real(w_0), imag(w_0)))
plot(x, real(phi_1), label="%.3f %.3f" % (real(w_1), imag(w_1)))
plot(x, real(phi_2), label="%.3f %.3f" % (real(w_2), imag(w_2)))
gkcStyle.plotZeroLine(min(kx), max(kx))
xlim((-5, 5.))

xlabel("$X$")
ylabel("$\phi(x)$")

####
subplot(122)
plot(x, imag(phi_0), label="%.3f %.3f" % (real(w_0), imag(w_0)))
plot(x, imag(phi_1), label="%.3f %.3f" % (real(w_1), imag(w_1)))
plot(x, imag(phi_2), label="%.3f %.3f" % (real(w_2), imag(w_2)))
gkcStyle.plotZeroLine(min(kx), max(kx))

xlim((-5, 5.))

xlabel("$X$")
ylabel("$\phi(x)$")

legend(loc='best').draw_frame(0)

savefig("ModeStructure.pdf", bbox_inches='tight')



