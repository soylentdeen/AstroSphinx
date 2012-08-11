.. _guide:

Uncle!
======

I have been attempting to implement the linear DELO (Diagonal Element Lambda Operator) algorithm of Rees et al. (1989). The algorithm seems to be relatively straight forward, but I have been unable to get it to converge fo the results of the Runge-Kutta algorithm.  I find the Runge-Kutta algorithm results reliable because it reproduces the reults of the scalar version of Scalar Moog to an accuracy of :math:`<0.1\%`.

The Runge-Kutta algorithm uses a canned FORTRAN 77 integration routine called dop853.f from `Ernst Hairer <http://www.unige.ch/~hairer/>`_, a well tested 8th order Runge-Kutta code.  The RK version of MoogStokes can be described in Pseudocode as::

    Read in Model Atmosphere
    Read in Line List
    Solve Boltzmann and Saha equations for each species
    Calculate line center opacities for each line/layer
    current wl = wl start
    while current wl < wl stop do
        for i = 0, ntau
            compute continuuous absorption quantities for layer i
            compute line absorption quatntities for layer i
            compute opacity matrix K for layer i
            compute emission vector J for layer i
        endfor
        initial Stokes vector = (B,0,0,0)^T
        initial Continuum = B

        Call RungeKutta routine dop853(Stokes, Continuum)

        I= Stokes(1)/Continuum
        Q= Stokes(2)/Continuum
        U= Stokes(3)/Continuum
        V= Stokes(4)/Continuum

        current wl += delta wl
    enddo

    Derivs(tau,I, Continuum) 
              ! <-Subroutine required by dop853.f to compute the 
              ! derivatives of each Stokes vector at a given tau
       dI = K(1,1,tau)*I(1) + K(1,2,tau)*I(2) + K(1,3,tau)*I(3) +
            K(1,4,tau)*I(4) - K(1,1,tau)*J(1,tau)
       dQ = K(2,1,tau)*I(1) + K(2,2,tau)*I(2) + K(2,3,tau)*I(3) +
            K(2,4,tau)*I(4) - K(2,1,tau)*J(2,tau)
       dU = K(3,1,tau)*I(1) + K(3,2,tau)*I(2) + K(3,3,tau)*I(3) +
            K(3,4,tau)*I(4) - K(3,1,tau)*J(3,tau)
       dV = K(4,1,tau)*I(1) + K(4,2,tau)*I(2) + K(4,3,tau)*I(3) +
            K(4,4,tau)*I(4) - K(4,1,tau)*J(4,tau)
       dContinuum = kaplam/kapref * (Continuum - J(1, tau)

       return(dI, dQ, dU, dV, dContinuum)

Because the model atmosphere is only defined on a rather coarse grid (:math:`\delta \tau \sim 0.1` dex), and the RK algorithm often samples the system of coupled differential equations at much smaller steps (:math:`\delta \tau \sim 0.001` dex), the DERIVS function must interpolate between the layers defined in the model atmosphere.

Here is the source code which calculates the opacity matrix :math:`\vec{\kappa}` as a function of optical depth through the photosphere for a single wavelength::

  FORTRAN CODE
         do i=1,ntau
            phi_I=(phi_opacity(i,2)*sin(phi_angle)**2.0+
        .                (phi_opacity(i,1)+phi_opacity(i,3))*(1.0+
        .                cos(phi_angle)**2.0)/2.0)/2.0
            phi_Q=(phi_opacity(i,2)-(phi_opacity(i,1)+phi_opacity(i,3))
        .                /2.0)*sin(phi_angle)**2.0*cos(2.0*chi_angle)/2.0
            phi_U=(phi_opacity(i,2)-(phi_opacity(i,1)+phi_opacity(i,3))
        .                /2.0)*sin(phi_angle)**2.0*sin(2.0*chi_angle)/2.0
            phi_V=(phi_opacity(i,1)-phi_opacity(i,3))*cos(phi_angle)/2.0
            psi_Q=(psi_opacity(i,2)-(psi_opacity(i,1)+psi_opacity(i,3))
        .                /2.0)*sin(phi_angle)**2.0*cos(2.0*chi_angle)/2.0
            psi_U=(psi_opacity(i,2)-(psi_opacity(i,1)+psi_opacity(i,3))
        .                /2.0)*sin(phi_angle)**2.0*sin(2.0*chi_angle)/2.0
            psi_V=(psi_opacity(i,1)-psi_opacity(i,3))*cos(phi_angle)/2.0

   c*****  The total opacity (line+continuum)
            kaptot(i) = kaplam(i) + phi_I
   c*****  Assemble the Opacity matrix (K)
            kappa(1,1,i)=kaptot(i)/kapref(i)
            kappa(1,2,i)=phi_Q/kapref(i)
            kappa(1,3,i)=phi_U/kapref(i)
            kappa(1,4,i)=phi_V/kapref(i)
            kappa(2,1,i)=phi_Q/kapref(i)
            kappa(2,2,i)=kaptot(i)/kapref(i)
            kappa(2,3,i)=psi_V/kapref(i)
            kappa(2,4,i)=(-1.0*psi_U)/kapref(i)
            kappa(3,1,i)=phi_U/kapref(i)
            kappa(3,2,i)=(-1.0*psi_V)/kapref(i)
            kappa(3,3,i)=kaptot(i)/kapref(i)
            kappa(3,4,i)=psi_Q/kapref(i)
            kappa(4,1,i)=phi_V/kapref(i)
            kappa(4,2,i)=psi_U/kapref(i)
            kappa(4,3,i)=(-1.0*psi_Q)/kapref(i)
            kappa(4,4,i)=kaptot(i)/kapref(i)

   c*****  Assumes LTE for the Source Function
            source = Planck(t(i))

   c*****  Assembles the Emission matrix (J')
            emission(1,i)=source!*kaptot(i)/kapref(i)
            emission(2,i)=source!*phi_Q/kapref(i)
            emission(3,i)=source!*phi_U/kapref(i)
            emission(4,i)=source!*phi_V/kapref(i)
   
            eta0(i) = kaplam(i)/kapref(i)
         enddo


I originally attempted a linear interpolation scheme between the neighboring boundary layers in the model atmosphere.  This seemed to do OK, reproducing the results of the Scalar Moog to :math:`\sim 0.01\%`.

I also attempted cubic spline interpolation.  Deep in the atmosphere, the I/C ratio oscillates a little bit around the linearly interpolated ratio, but it quickly converges to the same value.

DELO Algorithm
==============

The Runge-Kutta algorithm is an accurate way to compute the emergent Stokes vectors.  However, the RK algorithm is also processor intensive, and is not suitable for the calculation of a grid dense enough for my purposes.  For this reason, I investigate the DELO (Diagonal Element Lambda Operator) algorithm described by Rees et al. (1989).  The DELO algorithm is significanly faster than the Runge-Kutta algorithm with a modest tradeoff in accuracy.

The catch is that I have yet to be able to coax the DELO algorithm into converging to the same results as the Runge-Kutta routine, no matter how fine I make the regridding.::

    Read in Model Atmosphere
    Read in Line List
    Solve Boltzmann and Saha equations for each species
    Calculate line center opacities for each line/layer
    current wl = wl start
    while current wl < wl stop do
        for i = 0, ntau
            compute continuuous absorption quantities for layer i
            compute line absorption quatntities for layer i
            compute opacity matrix K for layer i
            compute emission vector J for layer i
        endfor
        initial Stokes vector = (B,0,0,0)^T
        initial Continuum = B

        Call DELO routine dop853(Stokes, Continuum)

        I= Stokes(1)/Continuum
        Q= Stokes(2)/Continuum
        U= Stokes(3)/Continuum
        V= Stokes(4)/Continuum

        current wl += delta wl
    enddo

the DELO algorithm is described in detail in Rees et al. (1989).  I shall repeat it here:

The opacity matrix :math:`\vec{\kappa}` is given by 

.. math::
   \vec{\kappa} = \kappa_C \vec{1} + \kappa_0 \vec{\phi}


.. math::
   \vec{\kappa} = \left(\begin{array}{cccc}
                   \kappa_C & 0 & 0 & 0 \\
                   0 & \kappa_C & 0 & 0 \\
                   0 & 0 & \kappa_C & 0 \\
                   0 & 0 & 0 & \kappa_C \end{array} \right)
           + \kappa_0\left(\begin{array}{cccc}
             \phi_I & \phi_Q & \phi_U & \phi_V \\
             \phi_Q & \phi_I & \psi_V & \psi_U \\
             \phi_U & -\psi_V & \phi_I & \psi_Q \\
             \phi_V & \psi_U & -\psi_Q & \phi_I \end{array} \right)

.. math::
   \vec{\kappa} = \left(\begin{array}{cccc}
            \kappa_I & \kappa_0\phi_Q & \kappa_0\phi_U & \kappa_0\phi_V \\
            \kappa_0\phi_Q & \kappa_I & \kappa_0\psi_V & \kappa_0\psi_U \\
            \kappa_0\phi_U & -\kappa_0\psi_V & \kappa_I & \kappa_0\psi_Q \\
            \kappa_0\phi_V & \kappa_0\psi_U & -\kappa_0\psi_Q & \kappa_I \\
            \end{array} \right)

where :math:`\kappa_I = \kappa_C + \kappa_0\phi_0`.  The DELO algorithm uses a modified absorption matrix :math:`\vec{\kappa'} = \frac{\vec{\kappa}}{\kappa_I} - \vec{1}`, which results in a matrix with zeroes along its diagonal.

.. math::
   \vec{\kappa'} = \left(\begin{array}{cccc}
    0 & \frac{\kappa_0\phi_Q}{\kappa_I} & \frac{\kappa_0\phi_U}{\kappa_I} & \frac{\kappa_0\phi_V}{\kappa_I} \\
    \frac{\kappa_0\phi_Q}{\kappa_I} & 0 & \frac{\kappa_0\psi_V}{\kappa_I} & \frac{\kappa_0\psi_U}{\kappa_I} \\
    \frac{\kappa_0\phi_U}{\kappa_I} & \frac{-\kappa_0\psi_V}{\kappa_I} & 0 & \frac{\kappa_0\psi_Q}{\kappa_I} \\
    \frac{\kappa_0\phi_V}{\kappa_I} & \frac{\kappa_0\psi_U}{\kappa_I} & \frac{-\kappa_0\psi_Q}{\kappa_I} & 0 \\
            \end{array} \right)

The ordinary source function :math:`\vec{J}` is

.. math::
   \vec{J} =\kappa_C S_C \left(\begin{array}{c} 1 \\ 0 \\ 0 \\ 0 \\ \end{array} \right) + \kappa_0 S_L \left(\begin{array}{cccc}
             \phi_I & \phi_Q & \phi_U & \phi_V \\
             \phi_Q & \phi_I & \psi_V & \psi_U \\
             \phi_U & -\psi_V & \phi_I & \psi_Q \\
             \phi_V & \psi_U & -\psi_Q & \phi_I \end{array} \right) +
        \left(\begin{array}{c} 1 \\ 0 \\ 0 \\ 0 \\ \end{array} \right)

where in the assumption of LTE, the continuum and line source functions are equal to the Planck Function :math:`\left(S_C = S_L = B\right)`

The DELO algorithm modifies he source function by dividing it by :math:`\kappa_I` as well, giving

.. math::
   \vec{S'} = \frac{\vec{J}}{\kappa_I} = B\left(\begin{array}{c} 1 \\
            \frac{\kappa_0\phi_Q}{\kappa_I} \\
            \frac{\kappa_0\phi_U}{\kappa_U} \\
            \frac{\kappa_0\phi_V}{\kappa_V} \\ \end{array} \right)

Following the prescription of Rees et al. (1989), the transfer euqation can be written as:

.. math::
   \frac{d\vec{I}}{d\tau} = \vec{I} - \left(\vec{S}' - \vec{K}'\vec{I}\right)

This in turn can be solved by a lambda transform, and after some algebra, a linear equation appears:

.. math::
   \vec{X_i}\cdot\vec{I}\left(\tau_i\right) = \vec{Y_i}\cdot\vec{I}\left(\tau_{i+1}\right)+\vec{Z_i}

where

.. math::
   \vec{X_i} &= \vec{1} + \left(\alpha_i - \beta_i \right)\vec{K}'_i \\
   \vec{Y_i} &= \left(\epsilon_i \vec{1} - \beta_i\vec{K}'_{i+1}\right) \\
   \vec{Z_i} &= \left(\alpha_i - \beta_i\right)\vec{S}'_i + \beta_i\vec{S}'_{i+1} \\
   \alpha_i &= 1 - \epsilon_i \\
   \beta_i &= \frac{\left(1-\left(1+\delta_i\right)\epsilon_i\right)}{\delta_i} \\
   \epsilon_i &= e^{-\delta_i} \\
   \delta_i &= \tau_{i+1} - \tau_i \\

While the equations look complicated, they are actually rather simple, and are easily translated into FORTRAN77 code::

  DELO Routine:
        Stokes(1) = emission(1,ntau)
        Stokes(2) = dble(0.0)
        Stokes(3) = dble(0.0)
        Stokes(4) = dble(0.0)
        continuum = Stokes(1)

  c**** Fill in initial emission, tau, and kappa values before loop starts
  
        call dcopy(4, emission(:,ntau), 1, emiss_interp(:,1), 1)
        tau_interp(1) = tauref(ntau)*kaptot(ntau)/kapref(ntau)
        tau_interp_c(1) = tauref(ntau)*kaplam(ntau)/kapref(ntau)
        call interp_opacities(log10(tauref(ntau)),
       .        kappa_interp, 1, emiss_interp, 1, tau_interp,tau_interp_c)
        kappa_order(1) = 1
        kappa_order(2) = 2
        emiss_order(1) = 1
        emiss_order(2) = 2

  c**** set tau stepsize
        delta_tau = -0.01
  
        do logtau=log10(tauref(ntau))+delta_tau,
       .              log10(tauref(1)),delta_tau
  c****      Get interpolated emission, tau, and kappa values for this step

           call interp_opacities(logtau, kappa_interp,
       .        kappa_order(2), emiss_interp, emiss_order(2), tau_interp,
       .        tau_interp_c)

  c****     Calculates change in total optical depth (continuum+line)

           dtau = (tau_interp(emiss_order(1))-tau_interp(emiss_order(2)))
       .              *cos(viewing_angle)
           etau = 2.71828183**(-dtau)

  c****      Calculates quantities for and constructs matrices

           alph = 1.0-etau
           bet =(1.0-(1.0+dtau)*etau)/dtau

           call dcopy(16,ones, 1, matX, 1)
           call dcopy(16,ones, 1, matY, 1)
           call daxpy(16,dble(alph-bet),kappa_interp(:,:,kappa_order(2)),
       .              1,matX,1)
           call dscal(16,etau, matY,1)
           call daxpy(16,dble(-1.0*bet),kappa_interp(:,:,kappa_order(1)),
       .              1,matY,1)

           call dcopy(4, emiss_interp(:,emiss_order(2)), 1, matS1, 1)
           call dcopy(4, emiss_interp(:,emiss_order(1)), 1, matZ, 1)
           call dscal(4, alph-bet, matS1, 1)
           call dscal(4, bet, matZ, 1)
           call daxpy(4, dble(1.0), matS1, 1, matZ, 1)

  c****      calculate the RHS of the equation.  Store in matZ
           call dgemv('N',4,4,dble(1.0),matY,4,Stokes,1,dble(1.0),matZ,1)

  c****      Solve the system of differential equations
           call dgesv(4,1,matX,4,IPIV,matZ,4,INFO)
  
           call dcopy(4, matZ, 1, Stokes, 1)
  
  c****     Now do the same thing for the continuum
           dtau=(tau_interp_c(emiss_order(1))-
       .         tau_interp_c(emiss_order(2)))*cos(viewing_angle)
           etau = 2.71828183**(-dtau)
           alph = 1.0 - etau
           bet =(1.0-(1.0+dtau)*etau)/dtau
           continuum=etau*continuum+(alph-bet)*
       .              emiss_interp(1,emiss_order(2))
       .             +bet*emiss_interp(1,emiss_order(1))
  
           if (kappa_order(1).eq.1)then
               kappa_order(1) = 2
               kappa_order(2) = 1
           else
               kappa_order(1) = 1
               kappa_order(2) = 2
           endif
           if (emiss_order(1).eq.1) then
               emiss_order(1) = 2
               emiss_order(2) = 1
           elseif (emiss_order(1).eq.2) then
               emiss_order(1) = 1
               emiss_order(2) = 2
           endif
        enddo
  c***   Once loop has ended, integration has reached edge of photosphere,
  c***   the emergent values are stored in Stokes(I,Q,U,V) and Continuum

I have checked, double-checked, triple-checked, and checked to the n-th degree that the logic contained in this algorithm.  I have ensure that I have been careful to interpolate all values correctly.  Still, I end up with the following comparison:



+---+
|   |
+---+


.. :Authors:
.. :Copyright:
