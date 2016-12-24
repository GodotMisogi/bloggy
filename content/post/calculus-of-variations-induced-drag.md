+++
date = "2016-12-24T12:30:48+05:30"
title = "Calculus of Variations - Induced Drag"
draft = false

+++

While reading through John D. Anderson Jr.'s derivation of minimum induced drag, I saw a cool application of the calculus of variations in one of the equations to deduce the required condition.

The equation that determines the downwash at a point is:

$$w(y\_0) = \frac{1}{4\pi}\int^{b/2}\_{-b/2} \frac{(\mathrm{d}\Gamma/\mathrm{d}y)}{y\_0 - y}\mathrm{d}y = \frac{1}{4\pi}\int^{b/2}\_{-b/2} \mathcal{L}(\Gamma,\Gamma',y)\;\mathrm{d}y$$ 

This effectively implies that the downwash can be expressed as a *functional* of $\Gamma$, i.e. $w\left[\Gamma(y)\right]$, and one can find the functional derivative to find the extremal point. The Euler-Lagrange equations thus take the following form:

$$ \frac{\partial{\mathcal{L}}}{\partial{\Gamma}} - \frac{\mathrm{d}}{\mathrm{d}y}\left(\frac{\partial{\mathcal{L}}}{\partial{\Gamma'}}\right) = 0 $$

Since there is no explicit $\Gamma$ dependence, the E-L equations simply become: $\partial{\mathcal{L}}/\partial{\Gamma'} = c$. This results in the following expression:

$$ y = y\_0 - \frac{1}{c} = k$$