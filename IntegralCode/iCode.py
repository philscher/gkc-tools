import numpy as np
import iCodeModule


class Species:
  def __init__(self, m=1., q=1., T=1., n=0., eta=0., name="Unamed"):
    self.m    = m
    self.q    = q
    self.T    = T
    self.n    = n
    self.eta  = eta
    self.name = name  


class solveDispersion:
  
  def __init__(self, N, verbose=True):

      # Mode / Number
      self.N       = N
      self.verbose = verbose
      
      # Setup Matrix
      self.L = np.zeros((N,N), dtype=complex)
    
 
  def solveEquation(self, w_init, ky, Setup, onlyValue=False, onlyDet=False, onlyWithQL=False):
        
        # Set basic parameters
        Ls, Ln, Lx, Debye2, N = Setup['Ls'], Setup['Ln'], Setup['Lx'], Setup['Debye2'], self.N
        dx, kx_list, dk  = Lx/N, 2.*np.pi/Lx * np.linspace(-(N-1)/2., (N-1)/2., N), 2.*np.pi/Lx
        X                = np.linspace(-Lx/2., Lx/2., N)
  
        def setupMatrix(w):

           # Resest Matrix & Recalculate (using Cython & C++)
           self.L[:,:] = 0.
           iCodeModule.setupMatrixPy(Setup['species'], w, ky, X, kx_list, Ls, Ln, self.N, self.L, dk*dx, Debye2)

        ################ Find Eigenvalue ##################3
       
        def getDeterminant(w):

           w = complex(w)
       
           # Setup new Matrix
           setupMatrix(w)

           # Get Minium absolute eigenvalue. 
           # 
           # The non-linear eigenvalue problem obeys 
           # det(L) = 0. This requirement is equal
           # to an eigenvalue 0 for the linear eigenvalue
           # problem det(L'- lambda I) = 0, which is identical
           # to the non-linear requirenement once lambda = 0.
           # 
           # [ This is less sensitive (and thus numerically less demanding) 
           # than caluculating directly the determinant
           eigvals  = np.linalg.eigvals(self.L)
           idx      = np.argmin(abs(eigvals))
           residual = eigvals[idx]

           if self.verbose :
              print ky, " w   :  %.8f+%.8f j" % (np.real(complex(w)), np.imag(complex(w))) ,  "  Determinant : %.2e " %  abs(residual)

           return residual	
   
        if onlyDet == True : return getDeterminant(w_init)

        try :
           omega = self.solveMuller(getDeterminant, w_init, tol =1.e-7, ftol=1.e-7, maxsteps=32)
        except:
           omega = float('nan') + 1.j * float('nan')

        ################ Find Eigenfunction ##################3
        if onlyValue == True: return omega

        # solution found for w0, get solution vector
        setupMatrix(omega) 

        # We found our eigenvalue omega, now we use the
        # inverse iteration to find the closest eigenvector
        # to the eigenvalue
        L__lambda_I = self.L - omega * np.eye(N)

        # Start with "random" inital eigenvector
        b = (1.+1.j) * np.ones(N)

        # Convergence is fast thus large iteration number not required
        # However, how to best check the error ? 
        # e.g. A dot phi
        # Tested also, Rayleigh-Quotient Iteration but less successfull.
        # R-Q Iteration only for hermitian matrices ?!
        residual = 1.e99
        
        for n in range(128):

            # Rescale b to unit vector 
            b = b/np.real(np.sqrt(np.vdot(b,b)))

            # Solve for next b
            b = np.linalg.solve(L__lambda_I, b)

            # calculate residul L(phi) = r
            r = np.dot(self.L,b)
            residual = np.real(np.sqrt(np.vdot(r,r)))

            if (residual < 1.e-10) : break
        
        print("I-I Residual : %.2e " % residual )
        
        # We have to transform to FFTW format, which has
        # form of [ k=0, k=1, ..., k = N/2, k = -(N/2-1), ..., k=-1 ] 
        K = np.append(np.append(b[N/2], b[N/2+1:]), b[:N/2])
        R = np.fft.ifft(K)
        R = np.append(np.append(R[N/2], R[N/2+1:]), R[:N/2])

        # Correct for phase
        R = R * np.exp(-1.j*np.arctan2(np.imag(sum(K)),np.real(sum(K))))

        #rescale to max(real(K)) = 1
        R = R / max(np.real(R))
       

        # Calculate quasi-linear heat-flux as
        #
        # Q_ql = sum_(k_x) omega_i / (k_x^2+k_y^2) |f(k_x)_ky|
        #
        # use rescaled value for calculations, is this correct ?
        b_resc = b / max(abs(b))
        if onlyWithQL == True: 
           QL = np.sum(np.imag(omega) * abs(b_resc) / (kx_list**2 + ky**2)) * dk
           print "Quasi-linear heat flux : ", QL
           return omega, QL

        return omega, kx_list, b, X, R

  """
    Muller's method
    http://en.wikipedia.org/wiki/Muller's_method

    Extraced from mpmath with minor modifications

  """
  def solveMuller(self, f, x0, tol =1.e-7, ftol=1.e-7, maxsteps=64):

     x0 , x1 , x2  = x0[0], x0[1], x0[2]
     fx0, fx1, fx2 = f(x0), f(x1), f(x2)

     for n in range(maxsteps):
   
        # calculate divided differences
        fx2x1 = (fx1 - fx2) / (x1 - x2)
        fx2x0 = (fx0 - fx2) / (x0 - x2)
        fx1x0 = (fx0 - fx1) / (x0 - x1)
   
        w       = fx2x1 + fx2x0 - fx1x0
        fx2x1x0 = (fx1x0 - fx2x1) / (x0 - x2)
   
        # update
        x0, fx0 = x1, fx1
        x1, fx1 = x2, fx2
   
        # denominator should be as large as possible => choose sign
        r = np.sqrt(w**2 - 4.*fx2*fx2x1x0)
        if abs(w - r) > abs(w + r): r = -r
   
        x2 -= 2.*fx2 / (w + r)
        fx2 = f(x2)
   
        error = abs(x2 - x1)
        if(error < tol and fx2**2 < ftol) : break

     # If root did not converges return nan
     if(error > tol or fx2**2 > ftol) : return float('nan')
   
     return x2
    

