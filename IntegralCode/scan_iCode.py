import sys
from iCode import *
from numpy import *

Nx  = 257
Nky = 128

#stdSetup_swITG = { 'Ln' : 1. , 'Ls' : 1./0.2, 'Lx' : 20., 'Debye2' : 0., \
# 'species' : [ Species(m=0.,n=0.,q=-1, eta=0., name="Adiab"),\
#               Species(n=1.,m=1.,q=1.,eta=5.,name="Ion") ,\
#               Species(m=1./1837.,n=1.,q=-1, eta=5., name="Electron") ]}


label = 'EtaScan_'
################################# Scan ########################
if len(sys.argv[:]) == 1:
  


 setup_list = []

 stdSetup_swITG_Adiab = { 'Ln' : 1. , 'Ls' : 1./0.2, 'Lx' : 20., 'Debye2' : 0., 'species' : [ Species(m=0.,n=1.,q=-1,name="Adiab"),  Species(n=1.,m=1.,q=1.,eta=5.,name="Ion") ] }

 setup_list.append(stdSetup_swITG_Adiab)
 
 for setup in setup_list:

  disp = solveDispersion(Nx)
  ky_list = logspace( -1., 1., Nky)

  results = []

  ################### scan for the ITG-sw mode ################
  # We trace the modes starting from the maximum growing point 
  # of the std-ITG mode
  
  def followLine(ky_split, w0_init):
  
    idx_l = where(ky_list <= ky_split)[0][::-1]
    idx_u = where(ky_list  > ky_split)[0]
    idx   = append( idx_u, idx_l)

    # Get inverse index god bless 
    # (http://stackoverflow.com/questions/11649577/how-to-invert-a-permutation-array-in-numpy) 
    idx_i = np.argsort(idx)

    omega_list = []
    qalin_list = []

    for n in range(Nky):

       if n == 0 : w0 = w0_init
       # Continue backwards using neighbouring eigenvalue
       if n == len(idx_u) : w0 = omega_list[0]

       w0, ql  = disp.solveEquation((w0, w0+0.001, w0+0.002j), ky_list[idx][n], setup, onlyWithQL=True)
       omega_list.append(w0)
       qalin_list.append(ql)
       print n , "/", Nky , "  ", ky_list[idx][n]
    
    print "-#####->",  array(qalin_list)[idx_i] 
    # return in standard order
    return array(omega_list)[idx_i], array(qalin_list)[idx_i] 
  
  ################################### Scan ##################
 
  #### Investigate std-modes
  ky_split_std = 0.5
  for w0_init in [ -.24 + 0.12j, -0.4 + 0.05j , - 0.8 - 0.01j, - 1.2 - 0.04j ]:
    omega, ql = followLine(ky_split_std, w0_init)
    results.append( (omega, ql) )
  
  #### Investigate sw-modes
  ky_split_sw = 3.
  # Very nice -0.4 + 0.02 
  for w0_init in [ -0.25 + 0.05j, -0.5 + 0.01j, -0.4 + 0.02 ]:
    omega, ql = followLine(ky_split_sw, w0_init)
    results.append( (omega, ql) )

  ######################### Data Ouput follows ############ 
  
  # Create two dimensional array
  col = len(results)
  row = len(ky_list)

  # time 2 due to imaginary part
  A = zeros((row, 3*col+1)) 
  print shape(A[:,0]), shape(ky_list)
  A[:,0] = array(ky_list)
  # We have offset due to ky_list column
  for n in range(len(results)) : print "--> ", array(results[n][1])
  for n in range(len(results)) : A[:,3*n+1], A[:,3*n+2], A[:,3*n+3] = real(array(results[n][0])), imag(array(results[n][0])), array(results[n][1])

  #stdSetup_swITG = { 'Ln' : 1. , 'Ls' : 1./0.2, 'Lx' : 20., 'Debye2' : 0., 'species' : [ Species(m=0.,n=1.,q=-1,name="Adiab"),  Species(n=1.,m=1.,q=1.,eta=5.,name="Ion") ] }
  id_str = "Shear_Ql_Nkx_%i_%.2f" % (Nx, 1./setup['Ls'])
  for  n in range(len(setup['species'])):
      id_str +=  "_Eta_%.2f_Mass_%.5f" %    ( setup['species'][n].eta,\
                                              setup['species'][n].m)
  savetxt(label + id_str + ".txt", A) 


