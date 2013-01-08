from pylab import *
import gkcStyle

from iCode import *


def solveMode(ky_list, w_init_list, Nx, Setup):

  print "Nx"
  disp = solveDispersion(Nx)
  
  for (ky, w_init) in zip(ky_list, w_init_list) : 
    w0, kx_list, K, X, R = disp.solveEquation(w_init, ky, Setup)
    
    clf()

    fig = gkcStyle.newFigure(ratio='1.41:1', basesize=9)
        
    fig.suptitle("$w0$ = %.4f %.4fi  $\pm$ %.2e %.2e i" % (real(w0), imag(w0), 0., 0.))
        
    
    ###################### Plot Fourier Modes ##########
    plot(kx_list, real(K), 'r.-', label="real")
    plot(kx_list, imag(K),  '.-', label="imag", color=gkcStyle.color_indigo)
        
    xlim((min(kx_list), max(kx_list))) 

    xlabel("$k_x$")
    ylabel("$\phi(k_x)$")
    legend(ncol=2).draw_frame(0)
        
    
    savefig(str(ky) + "_Plot1.pdf", bbox_inches='tight')
    ################### Plot real modes ########3
    clf()
    
    plot(X, real(R), 'r.-', label='real')
    plot(X, imag(R),  '.-', label='imag', color=gkcStyle.color_indigo)
    
    xlim((min(X), max(X)))
    
    xlabel("$x$")
    ylabel("$\phi(x)$")
    legend(ncol=2).draw_frame(0)
        
    
    savefig(str(ky) + "_Plot2.pdf", bbox_inches='tight')
    
    ################ Plot Contour ###############
    
    clf()
    
    Ly = 2.*pi/ky
    y = linspace(0., Ly, 256)
    
    KxKy = zeros((Nx, 65), dtype=complex)
    KxKy[:,1] = R
    
    XY = np.fft.irfft(KxKy, axis=1, n=256)
    
    xlabel("$x$")
    ylabel("$y$")
    
    contourf(X, y, XY.T, 20, vmin=-abs(XY).max(), vmax=abs(XY).max())
    colorbar()
    
    savefig(str(ky) + "_Plot3.pdf", bbox_inches='tight')
    
    # append and normalize
    #sol.append(np.fft.ifft(b/abs(b).max()))

    # Write out real mode structure
    savetxt("IntegralCode_ky_" + str(ky) + ".txt", array([X, real(K), imag(K)]).T ) 


def plotContours():

    ky = 0.5
    R = linspace(-0.5, 0.5, 8)
    I = linspace(-.3, 0.3, 8)
    V = zeros((len(R),len(I)), dtype=complex)

    fig = figure(figsize=(30,10))
   
    n = 0
    for r in range(len(R)):
     for i in range(len(I)):

      A = zeros((Nx,Nx), dtype=complex)
      iCode.setupMatrixPy(species, R[r]+1.j*I[i], ky, X, kx_list, Ls, Ln, Nx, A, dk*dx, lambda_D2)
      #val =  getMinEigenvalue(A)
      (sign, logdet) = np.linalg.slogdet(A)
      val = sign * logdet
      V[r,i] = val
      #print "x, y", R[r], I[i] , " r : ", val
      print n, "/", len(R) * len(I)
      n = n+1
  
    """
    #subplot(131)
    #norm = mpl.colors.Normalize(vmin = -1., vmax = 1.)
    #contourf(R,I,real(V), 100, vmin=-1., vmax=1., norm = norm)
    xlabel("Real")
    ylabel("Imag")
    cb = colorbar()
    cb.set_clim(vmin=-1, vmax=1)

    #subplot(132)
    #contourf(R,I,imag(V), 100, vmin=-1., vmax=1.)
    #norm = mpl.colors.Normalize(vmin = -1., vmax = 1.)
    contourf(R,I,imag(V), 100, vmin=-1., vmax=1., norm = norm)
    xlabel("Real")
    ylabel("Imag")
    cb = colorbar()
    cb.set_clim(vmin=-1, vmax=1)
    
    subplot(133)
    """
    pcolor(R,I,log10(abs(V)))
   
    xlabel("Real")
    ylabel("Imag")
    cb = colorbar()
    #cb.set_clim(vmin=0., vmax=1)
    
    
    #pcolor(R,I,imag(V))
    savefig("Contour.png")
    
    #print "(Integral)  Solution  is w : ",w0
    #print "(local)     Solution  is w : ",w_Local

    idx_l = where(ky_list < 2.55)

def plotEvolution():

  #ky_list = linspace(0.3, pi, 17)
  ky_list = logspace(-1, log10(pi), 3)

  results = []
 
  """
  for w0 in [ -0.4 + 0.05j , -0.2 + 0.01j, - 0.8 - 0.01j, - 1.2 - 0.04j ]:
    res = []
    for n in range(len(ky_list)):
     w0 = solveDispersion(ky_list[n], (w0, w0+0.01, w0+0.001j), doPlots=False)
     res.append(w0)

    results.append(res) 
  """ 
  w0 = -0.4 + 0.05j
  for w0_init in [ -0.25 + 0.05j, -0.5 + 0.01j, -0.4 + 0.02 ]:
    res = []
    idx_l = where(ky_list <  2.55)
    idx_u = where(ky_list >= 2.55)
    idx = append(idx_u, idx_l[0][::-1])
    ky_list_mod = ky_list[idx]
    for n in range(len(ky_list)):
     print " ky ---> " , ky_list[idx][n]
     if (abs(ky_list[idx][n] - 2.55) < 0.1) : w0 = w0_init
     w0 = solveDispersion(ky_list[idx][n], (w0, w0+0.01, w0+0.001j), doPlots=False)
     res.append(w0)
    res = array(res)
    results.append(res[idx]) 
  
  results = array(results)
   
  
  #fig = gkcStyle.newFigure(ratio='1.41:1', basesize=9)
  fig = gkcStyle.newFigure(ratio='1.41:1', basesize=18)
  subplot(121)
  for res in results:
      semilogx(ky_list, imag(res), '-', label='imag')
  gkcStyle.plotZeroLine(min(ky_list), max(ky_list))
  
  xlim((min(ky_list), max(ky_list)))
  xlabel("$k_y$")
  ylabel("Growthrates $\\gamma$")

  #savefig("Growthrate.pdf", bbox_inches='tight')
 
  subplot(122)
  #clf()
  for res in results:
      semilogx(ky_list, real(res), '-', label='imag')
  gkcStyle.plotZeroLine(min(ky_list), max(ky_list))
  
  xlim((min(ky_list), max(ky_list)))
  xlabel("$k_y$")
  ylabel("Frequency $\\gamma$")
  savefig("Frequency.pdf", bbox_inches='tight')

  # Save max tendency
  savetxt("Tendency.txt", array([ky_list, real(results), imag(results)]).T)


############################## Settings for Integral Mode ######################################
# My setup
#species = [  Species(name= "Adiab"),   Species(1.,1.,1.,1.,5., "Ion"),  Species(1./1837.,-1.,1.,1., 2., "Electron") ]

stdSetup_swITG = { 'Ln' : 1. , 'Ls' : 1./0.2, 'Lx' : 20., 'Debye2' : 0., 'species' : [ Species(m=0.,n=1.,q=-1,name="Adiab"),  Species(n=1.,m=1.,q=1.,eta=5.,name="Ion") ] }

## Gao Setup
#species = [  Species(name= "Adiab"),   Species(m=1836.,q=1.,T=1.,n=1.,eta=0., name="Ion"),  Species(m=1.,q=-1.,T=1.,n=1., eta=3., name="Electron") ]
#Ln, Ls, Lx, Ly, lambda_D2, ky_list, w0 = 1., 0.025, 30., 256., 1., [0.3], -0.01 + 0.025j

w0 = -0.25 + 0.05j
solveMode([0.491], [(w0, w0-0.01, w0+0.001j)], 513,  stdSetup_swITG )
#plotEvolution()
#plotContours()


######################## Setup Grid ######################


