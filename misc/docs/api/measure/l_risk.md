# ℓ-risk (loss-based Bayes risk)

```{py:module} qif.measure.l_risk
```

`qif.measure.l_risk`

A loss function $\ell(w,x)$ gives the adversary's loss when she takes action $w \in \mathcal W$ and the secret is $x$.
The $\ell$-risk $U_\ell(\pi) = \min_w \sum_x \pi_x\, \ell(w,x)$ is her expected loss when she takes the best
action; it is an uncertainty measure, so leakage is computed as prior minus (or divided by) posterior risk.

As for {doc}`g_vuln`, a loss function is given either as a loss matrix `L`, with one row per action and one column per
secret, or as a function `l(w, x)`, in which case the actions are assumed to be the secrets.

## Functions

```{eval-rst}
.. autofunction:: prior

   The :math:`\ell`-risk :math:`U_\ell(\pi) = \min_w \sum_x \pi_x\, \ell(w,x)`.
```

```{eval-rst}
.. autofunction:: posterior

   The posterior :math:`\ell`-risk
   :math:`U_\ell[\pi\triangleright C] = \sum_y \min_w \sum_x \pi_x C_{x,y}\, \ell(w,x)`.
```

```{eval-rst}
.. autofunction:: add_leakage

   The additive leakage :math:`U_\ell(\pi) - U_\ell[\pi\triangleright C]`.
```

```{eval-rst}
.. autofunction:: mult_leakage

   The multiplicative leakage :math:`U_\ell(\pi) / U_\ell[\pi\triangleright C]`.
```

```{eval-rst}
.. autofunction:: strategy

   The optimal action for each observation: an array whose element :math:`y` is the :math:`w` minimizing
   :math:`\sum_x \pi_x C_{x,y}\, \ell(w,x)`.
```

```{eval-rst}
.. autofunction:: add_capacity

   The additive capacity for a fixed prior, over a class of loss functions. Since the additive leakage of :math:`\ell`
   is the same as that of the gain function :math:`k - \ell` (for any constant :math:`k`), this is equal to
   :py:func:`qif.measure.g_vuln.add_capacity`: the class contains the loss functions with
   :math:`\max U_\ell - \min U_\ell \le 1`, or, if ``one_spanning_g`` is true, those with
   :math:`\max \ell - \min \ell \le 1`.
```

```{eval-rst}
.. autofunction:: loss_to_gain

   Converts a loss function into a gain function :math:`g(w,x) = k - \ell(w,x)`, where :math:`k` is the maximum of
   :math:`\ell` over ``n_guesses`` actions and ``n_secrets`` secrets. Then :math:`V_g = k - U_\ell`.
```

```{eval-rst}
.. autofunction:: l_zero_one

   The 0-1 loss function: :math:`\ell(w,x) = 0` if :math:`w = x`, :math:`1` otherwise. :math:`U_\ell` is Bayes risk.
```

```{eval-rst}
.. autofunction:: L_zero_one

   The matrix of the 0-1 loss function on ``n_rows`` secrets.
```
