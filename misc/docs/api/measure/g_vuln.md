# g-vulnerability

```{py:module} qif.measure.g_vuln
```

`qif.measure.g_vuln`

A gain function $g(w,x)$ gives the adversary's gain when she takes action $w \in \mathcal W$ and the secret is $x$.
The $g$-vulnerability $V_g(\pi) = \max_w \sum_x \pi_x\, g(w,x)$ is her expected gain when she takes the best
action.

A gain function is given either as a gain matrix `G`, with one row per action and one column per secret
($G_{w,x} = g(w,x)$), or as a function `g(w, x)`, in which case the actions are assumed to be the secrets
($\mathcal W = \mathcal X$).

## Functions

```{eval-rst}
.. autofunction:: prior

   The :math:`g`-vulnerability :math:`V_g(\pi) = \max_w \sum_x \pi_x\, g(w,x)`.
```

```{eval-rst}
.. autofunction:: posterior

   The posterior :math:`g`-vulnerability
   :math:`V_g[\pi\triangleright C] = \sum_y \max_w \sum_x \pi_x C_{x,y}\, g(w,x)`.
```

```{eval-rst}
.. autofunction:: add_leakage

   The additive :math:`g`-leakage :math:`\mathcal L^+_g(\pi, C) = V_g[\pi\triangleright C] - V_g(\pi)`.
```

```{eval-rst}
.. autofunction:: mult_leakage

   The multiplicative :math:`g`-leakage :math:`\mathcal L^\times_g(\pi, C) = V_g[\pi\triangleright C] / V_g(\pi)`.
```

```{eval-rst}
.. autofunction:: strategy

   The optimal action for each observation: an array whose element :math:`y` is the :math:`w` maximizing
   :math:`\sum_x \pi_x C_{x,y}\, g(w,x)`.
```

```{eval-rst}
.. autofunction:: add_capacity

   The additive capacity for a fixed prior: the maximum of :math:`\mathcal L^+_g(\pi, C)` over all gain functions
   :math:`g` in a class.

   By default the class contains the one-bounded gain functions (:math:`g \le 1`), and the capacity is
   :math:`1 - \sum_y \min_{x \in \lceil\pi\rceil} C_{x,y}`, where :math:`\lceil\pi\rceil` is the support of ``pi``.
   If ``one_spanning_g`` is true, the class contains the gain functions whose range has size at most :math:`1`
   (:math:`\max g - \min g \le 1`).
```

```{eval-rst}
.. autofunction:: g_id

   The identity gain function :math:`g_\mathrm{id}(w,x) = 1` if :math:`w = x`, :math:`0` otherwise.
   :math:`V_{g_\mathrm{id}}` is Bayes vulnerability.
```

```{eval-rst}
.. autofunction:: G_id

   The matrix of :math:`g_\mathrm{id}` on ``n_rows`` secrets (the identity matrix).
```

```{eval-rst}
.. autofunction:: g_add

   A gain matrix :math:`G` such that :math:`V_G = V_{G_1} + V_{G_2}`. Its actions are the pairs of actions of
   ``G1`` and ``G2``.
```

```{eval-rst}
.. autofunction:: g_from_posterior

   A gain matrix :math:`G'` such that :math:`V_{G'}(\pi) = V_G[\pi\triangleright C]` for all priors :math:`\pi`. Its
   actions are the strategies :math:`S : \mathcal Y \to \mathcal W`, with :math:`G'_{S,x} = \sum_y C_{x,y}\, G_{S(y),x}`,
   so its size grows as :math:`|\mathcal W|^{|\mathcal Y|}`.
```

```{eval-rst}
.. autofunction:: g_to_bayes

   Expresses :math:`g`-vulnerability on the prior ``pi`` as Bayes vulnerability. Returns ``(rho, R, a, b)``, where
   ``rho`` is a distribution on actions and ``R`` a channel from actions to secrets, such that
   :math:`V_g[\pi\triangleright C] = a\, V_1[\rho\triangleright RC] + b` for every channel :math:`C`.
```
