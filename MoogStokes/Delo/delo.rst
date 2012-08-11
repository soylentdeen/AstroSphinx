.. _guide:

Polarized Radiative Transfer Through a Magentized Stellar Atmosphere
====================================================================

Radiative Transfer Through a Stellar Atmosphere
-----------------------------------------------

The emergent light from a star is a composite of all the outward-bound photons
which manage to escape from their last scattering without further interaction
with atoms in the stellar atmosphere.  Photons created or scattered deep within
the star have almost no chance of escape, while photons created or scattered
at the edge of the photosphere are virtually guaranteed to escape.  The
conditions within a stellar photosphere range from extremely hot and dense (at
the base of the photosphere) to cooler and extrememly tenuous (at the edge of
the photosphere).  Applying the equation of radiative transfer to the case of a
stellar photosphere:

.. math::
   \frac{d\vec{I}}{d\tau}=-\kappa\vec{I} + \vec{J}

where :math:`\vec{I}` is the intensity vector of the radiation, :math:`\tau` is
the optical depth, :math:`\kappa` is the opacity matrix, and :math:`\vec{J}` is
the emission vector.

Scalar Radiative Transfer
-------------------------

Moog (Sneden 1976) uses contribution functions to calculate the emergent
spectrum.  This is not possible for polarized radiation (at least not in this
incarnation).

Polarized Radiative Transfer
----------------------------

Unno (1956)
^^^^^^^^^^^
Unno (1956) first attacked the problem of radiative transfer through a magnetic 
medium, calculating the propogation of the Stokes vectors through a magnetized
atmosphere.  The Unno equations accounted for the absorption profiles due to the
:math:`\sigma` and :math:`\pi` components and generalized the equation of
radiative transfer to the four Stokes vectors.

Rachovsky (1962)
^^^^^^^^^^^^^^^^
Rachovsky (1962) slightly modified the equations of Unno to include the
magneto-optical effects (so-called anomalous dispersion) due to changes in the
imaginary component of the index of refraction.

Degl'Innocenti (1976)
^^^^^^^^^^^^^^^^^^^^^
E. Landi Degl'Innocenti (1976) developed a fourth-order Runge-Kutta integration
code (MALIP) to calculate the Stokes parameters profiles of magnetoactive
Fraunhofer lines.  The derivatives of each of the Stokes parameters are given by
the following equations:

.. math::
  \mu\frac{dI}{d\tau}&=\eta_0(I-B_T)+\eta_I(I-S)+\eta_QQ+\eta_U+\eta_VV\\
  \mu\frac{dQ}{d\tau}&=\eta_Q(I-S)+(\eta_0+\eta_I)Q+\rho_VU-\rho_UV\\
  \mu\frac{dU}{d\tau}&=\eta_U(I-S)-\rho_VQ+(\eta_0+\eta_I)U+\rho_QV\\
  \mu\frac{dV}{d\tau}&=\eta_V(I-S)+\rho_VQ-\rho_QU+(\eta_0+\eta_I)V\\

where:
* :math:`\mu=\cos \gamma`, (where :math:`\gamma` is the angle between the
  line of sight and the normal to the stellar surface)
* :math:`\tau` is the optical depth, counted at a particular reference
  wavelength (:math:`\lambda_{ref}=500` nm) as a function of which
  the model-atmosphere and the other physical quantities are specified
* :math:`\eta_0=\frac{k_{\lambda_0}^{(c)}}{k_{\lambda_ref}^{(c)}}` is the ratio
  between the continuum opacities at the current wavelength :math:`\lambda_0`
  and the reference wavelength of the model atmosphere :math:`\lambda_{ref}`.


The full opacity matrix :math:`\kappa` is given by the following:

.. math::
  \kappa &= \kappa_c\vec{1} + \kappa_0\vec{\Phi}

where:

.. math::
  \vec{\Phi} &= \left( \begin{array}{cccc}
                 \phi_I & \phi_Q & \phi_U & \phi_V \\
                 \phi_Q & \phi_I & \phi'_V & \phi'_U \\
                 \phi_U & -\phi'_V & \phi'_I & \phi_Q \\
                 \phi_V & \phi'_U & -\phi'_Q & \phi_I \end{array} \right)

and :math:`\kappa_c` is the continuum opacity, and :math:`\kappa_0` is the line
center opacity for the spectral line in question.  

Diagonal Element Lambda Operator
--------------------------------

In order to trace the Stokes parameters through the stellar atmosphere, I use
the Diagonal Element Lambda Operator method developed by Rees & Murphy (1989).

Described in Rees & Murphy (1989), the DELO method works by 


Current State of Affairs
------------------------

I have an algorithm which works using a Runge-Kutta integrator.  This produces
the correct value for the non-magnetic case (as compared to the scalar Moog 
code) don't know if it produces the correct answer when the magnetic field is
turned on.  This will require me to synthesize an entire stellar disk.
Synthesizing an entire grid with this algorithm will be time consuming to say
the least.

I have attempted to implement the DELO algorithm as described in 

+---+
|   |
+---+

.. :Authors:
.. :Copyright:
