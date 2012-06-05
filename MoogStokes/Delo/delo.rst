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

Unno (1956) first attacked the problem of radiative transfer through a magnetic 
medium, calculating the propogation of the Stokes vectors through a magnetized
atmosphere.  The Unno equations accounted for the absorption profiles due to the
:math:`\sigma` and :math:`\pi` components and generalized the equation of
radiative transfer to the four Stokes vectors.

Rachovsky (1962) slightly modified the equations of Unno to include the
magneto-optical effects (so-called anomalous dispersion) due to changes in the
imaginary component of the index of refraction.  The full opacity matrix
:math:`\kappa` is given by the following:

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


+---+
|   |
+---+

.. :Authors:
.. :Copyright: