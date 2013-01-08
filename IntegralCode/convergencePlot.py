from pylab import *
import gkcStyle

fig = gkcStyle.newFigure(ratio="2.33:1", basesize=12)

Nx = array([ 33., 65., 127., 197., 257., 369., 513., 713.])
omega = [ -0.260 - 0.004j , - 0.218 + 0.061j, -0.229 + 0.082j, -0.233 + 0.0882j, - 0.235 + 0.0905j, -0.237 + 0.09267j, -0.238 + 0.094j, -0.239 + 0.0949j]

# Frequency
subplot(211)
plot(1./Nx, real(omega), 'o-')

xticks([0., 0.005, 0.01, 0.02, 0.04], [" \n$\\infty$", "200", "100", "50", "25"])

xlabel("Nx")
ylabel("Frequency")

subplot(212)

plot(1./Nx, imag(omega), 'o-')
xticks([0., 0.005, 0.01, 0.02, 0.04], [" \n$\\infty$", "200", "100", "50", "25"])

xlabel("Nx")
ylabel("Frequency")


savefig("Convergence.pdf", bbox_inches='tight')


