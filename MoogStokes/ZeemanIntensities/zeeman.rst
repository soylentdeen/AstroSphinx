.. _guide:

Zeeman Intensities
==================

Calculation of Zeeman Intensities
---------------------------------

Condon and Shortley (1967) give the following equations for calculating the relative intensities of Zeeman components:

.. math::
   <\alpha, J, M \| T \| \alpha, J+1, M\pm 1>^2 &= \frac{1}{4} \left(J\pm M+1\right)\left(J\pm M+2\right) \left(\bold{\hat{i}}\pm \bold{\hat{j}} \right)\\
   <\alpha, J, M \| T \| \alpha, J+1, M>^2 &= \left(J+1\right)^2 - M^2 \bold{\hat{k}}\\
   <\alpha, J, M \| T \| \alpha, J, M\pm 1>^2 &= \frac{1}{4} \left(J\mp M\right)\left(J\pm M+1\right) \left(\bold{\hat{i}}\pm \bold{\hat{j}} \right)\\
   <\alpha, J, M \| T \| \alpha, J, M>^2 &= M^2 \bold{\hat{k}}\\
   <\alpha, J, M \| T \| \alpha, J-1, M\pm 1>^2 &= \frac{1}{4} \left(J\mp M\right)\left(J\mp M-1\right) \left(\bold{\hat{i}}\pm \bold{\hat{j}} \right)\\
   <\alpha, J, M \| T \| \alpha, J-1, M>^2 &= J^2 - M^2 \bold{\hat{k}}

Case A
------
My naive assumption is to divide up the opacity contained in the original (unsplit) line between the Zeeman components by adjusting
their oscillator strengths accordingly, normalizing so that the sum of the oscillator strengths of all the Zeeman components equals
the oscillator strength of the unsplit line.

For example, consider the following spectral line: :math:`\lambda=11991.568\AA`, Si, :math:`\log gf=-0.109`.  For this transition, the
upper and lower states have total angular momentum numbers of :math:`J=0, J=1`, respectively.  The Zeeman effect will split the lower
state into three sublevels (:math:`m_j=-1,0,1`), while the upper state will remain a single level (:math:`m_j=0`).  In the simplest of all
configurations, the initial line splits into three Zeeman components (:math:`\Delta m_j=-1,0,1`).  Calculating the relative intensities of
the three components, the two sigma (:math:`\Delta m_j = \pm 1`) components should have a relative
intensity of :math:`\frac{1}{4}\left(0+0+1\right)\left(0+0+2\right)=\frac{1}{2}`.  The Pi component (:math:`\Delta m_j = 0`) should have a
relative intensity of :math:`\left(0+1\right)^2 - 0^2=1`.

I interpret this to mean that the :math:`\pi` component is twice as strong as either of the :math:`\sigma` components.  In order to 
reflect this, I adjust the oscillator strength of the :math:`\pi` component to be :math:`\log gf=-0.410`, and the oscillator strength 
of each of the two :math:`\sigma` components to be :math:`\log gf=-0.711`.  This way of dividing up the oscillator strengths ensures 
that the :math:`\pi` component is twice as strong as the :math:`\sigma` components, and that the total opacity in the three lines 
is the same as the opacity in the unsplit line.

However, when calculating the elements of the opacity matrix, this does not seem to asymptote to the correct value of the opacity at 
the limiting case of :math:`B=0`:

.. math::
   \Phi_I &= \frac{1}{2} \phi_p \sin^2 \gamma + \frac{1}{4}\left(\phi_r+\phi_b\right)\left(1+\cos^2\gamma\right)\\
   \Phi_Q &= \frac{1}{2} \left(\phi_p-\frac{1}{2}\left(\phi_r+\phi_b\right)\right)\sin^2\gamma\cos 2\chi\\
   \Phi_U &= \frac{1}{2} \left(\phi_p-\frac{1}{2}\left(\phi_r+\phi_b\right)\right)\sin^2\gamma\sin 2\chi\\
   \Phi_V &= \frac{1}{2} \left(\phi_r-\phi_b\right)\cos \gamma\\
   \Psi_Q &= \frac{1}{2} \left(\psi_p-\frac{1}{2}\left(\psi_r+\psi_b\right)\right)\sin^2\gamma\cos 2\chi\\
   \Psi_U &= \frac{1}{2} \left(\psi_p-\frac{1}{2}\left(\psi_r+\psi_b\right)\right)\sin^2\gamma\sin 2\chi\\
   \Psi_V &= \frac{1}{2} \left(\psi_r-\psi_b\right)\cos \gamma

where :math:`\phi_{p,b,r}` are the absorption profiles and :math:`\psi_{p,b,r}` are the anomalous dispersion profiles for the unshifted 
:math:`\pi` component, and blue/red-shifted :math:`\sigma` components, and :math:`\gamma, \chi` are angles related to the magnetic field
vector.

In the limiting case where :math:`B=0`, and :math:`\gamma=\chi=0`, we are looking at the center of a stellar disk of an unmagnetized star.
Since :math:`B=0`, the wavelength shift due to the Zeeman effect is zero, and all Zeeman components lie on top of each other, and the
spectrum should be identical to the spectrum calculated by a scalar code (i.e. Moog).

In Case A, the absorption and anomalous dispersion profiles differ for the different Zeeman components (i.e.
:math:`\phi_p=\frac{\phi_o}{2}, \phi_r=\phi_b=\frac{\phi_o}{4}, \psi_p=\frac{\psi_o}{2}, \psi_r=\psi_b=\frac{\psi_o}{4}`), so the elements
of the opacity matrix become:

.. math::
   \Phi_I &= \frac{1}{2} \frac{\phi_o}{2} \sin^2\left(0\right) + \frac{1}{4}\left(\frac{\phi_o}{4}+\frac{\phi_o}{4}\right)\left(1+\cos^2\left(0\right)\right)\\
          &= \frac{1}{4} \phi_o 0 + \frac{1}{4}\left(\frac{\phi_o}{2}\right)\left(1+1\right)\\
          &= \frac{\phi_o}{4}\\
   \Phi_Q &= \frac{1}{2} \left(\frac{\phi_o}{2}-\frac{1}{2}\left(\frac{\phi_o}{4}+\frac{\phi_o}{4}\right)\right)\sin^2\left(0\right)\cos\left(0\right)\\
          &= 0\\
   \Phi_U &= \frac{1}{2} \left(\frac{\phi_o}{2}-\frac{1}{2}\left(\frac{\phi_o}{4}+\frac{\phi_o}{4}\right)\right)\sin^2\left(0\right)\sin\left(0\right)\\
          &= 0\\
   \Phi_V &= \frac{1}{2}\left(\frac{\phi_o}{4}-\frac{\phi_o}{4}\right)\cos\left(0\right)\\
          &= 0\\
   \Psi_Q &= \frac{1}{2} \left(\frac{\psi_o}{2}-\frac{1}{2}\left(\frac{\psi_o}{4}+\frac{\psi_o}{4}\right)\right)\sin^2\left(0\right)\cos\left(0\right)\\
          &= 0\\
   \Psi_U &= \frac{1}{2} \left(\frac{\psi_o}{2}-\frac{1}{2}\left(\frac{\psi_o}{4}+\frac{\psi_o}{4}\right)\right)\sin^2\left(0\right)\sin\left(0\right)\\
          &= 0\\
   \Psi_V &= \frac{1}{2}\left(\frac{\psi_o}{4}-\frac{\psi_o}{4}\right)\cos\left(0\right)\\
          &= 0\\

All elements in the opacity matrix are zero, except for the :math:`\Phi_I=\frac{\phi_o}{4}`, only a fourth of what I would expect it to be.

Case B
------

However, Synthmag and other polarized radiative transfer codes don't seem to divide up the opacity this way.  Equation 11 in 
Rees et al. (1989) says that "the corresponding strength of this component is :math:`S_{i_j}` where :

.. math::
   \sum\limits_{i_j=1}^{N_j}S_{i_j} = 1

for :math:`j=r,p,b`, where :math:`r,b` correspond to the :math:`\pm \sigma` components, and :math:`p` correspond to the :math:`\pi` 
components.  If I'm interpreting this correctly, the normalization of oscillator strengths should be done three times (corresponding
to :math:`\Delta m_j = -1, 0, 1`).  In other words, the sum of oscillator strengths for all :math:`\pi` components should equal the 
original oscillator strength for the unsplit line.  The sum of oscillator strengths for all :math:`\sigma_+` components should also equal
the original oscillator strength for the unsplit line.  Similarly, the sum of oscillator strengths for all :math:`\sigma_-` components
should also equal the original oscillator strength for the unsplit line.  So, for the example described above, the oscillator strength for
each of the Zeeman components of the Si :math:`\lambda=11991.568\AA` line will be :math:`\log gf = -0.109`.

While I don't particularly understand the physical reason for this "normalization," I can at least see how it fits into the
calculations of the elements :math:`\left(\Phi_{I,Q,U,V},\Psi_{Q,U,V}\right)` of the opacity matrix:


In Case B, the absorption and anomalous dispersion profiles for each Zeeman component are the same (i.e.
:math:`\phi_p=\phi_r=\phi_b=\phi_o, \psi_p=\psi_r=\psi_b`=\psi_o`), so the elements of the opacity matrix become:

.. math::
   \Phi_I &= \frac{1}{2} \phi_o \sin^2\left(0\right) + \frac{1}{4}\left(\phi_o+\phi_o\right)\left(1+\cos^2\left(0\right)\right)\\
          &= \frac{1}{2} \phi_o 0 + \frac{1}{4}\left(2\phi_o\right)\left(1+1\right)\\
          &= \phi_o\\
   \Phi_Q &= \frac{1}{2} \left(\phi_o-\frac{1}{2}\left(\phi_o+\phi_o\right)\right)\sin^2\left(0\right)\cos\left(0\right)\\
          &= 0\\
   \Phi_U &= \frac{1}{2} \left(\phi_o-\frac{1}{2}\left(\phi_o+\phi_o\right)\right)\sin^2\left(0\right)\sin\left(0\right)\\
          &= 0\\
   \Phi_V &= \frac{1}{2}\left(\phi_o-\phi_o\right)\cos\left(0\right)\\
          &= 0\\
   \Psi_Q &= \frac{1}{2} \left(\psi_o-\frac{1}{2}\left(\psi_o+\psi_o\right)\right)\sin^2\left(0\right)\cos\left(0\right)\\
          &= 0\\
   \Psi_U &= \frac{1}{2} \left(\psi_o-\frac{1}{2}\left(\psi_o+\psi_o\right)\right)\sin^2\left(0\right)\sin\left(0\right)\\
          &= 0\\
   \Psi_V &= \frac{1}{2}\left(\psi_o-\psi_o\right)\cos\left(0\right)\\
          &= 0\\

So, in the limit of a non-magnetic star, the only non-zero element in the opacity matrix is :math:`\Phi_I=\phi_o`, which appears to approach
the correct value.




+---+
|   |
+---+

.. :Authors:
.. :Copyright:
