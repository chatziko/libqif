# Vulnerability and leakage measures

```{py:module} qif.measure
```

`qif.measure`

Each submodule implements a family of measures, with functions following the same naming scheme. For a prior
$\pi$ and a channel $C$ (see {doc}`../channel/index`):

`prior(pi)`
: The prior vulnerability $V(\pi)$ (or uncertainty $U(\pi)$) of the secret.

`posterior(pi, C)`
: The posterior vulnerability $V[\pi\triangleright C]$, that is the expected vulnerability of the posterior
  distributions, weighted by their probability: $\sum_y p(y)\, V(p_{X|y})$ (similarly for $U$).

`add_leakage(pi, C)`, `mult_leakage(pi, C)`
: The additive and multiplicative leakage, $V[\pi\triangleright C] - V(\pi)$ and $V[\pi\triangleright C] / V(\pi)$.
  For uncertainty measures the order is reversed, $U(\pi) - U[\pi\triangleright C]$ and
  $U(\pi) / U[\pi\triangleright C]$, so that leakage is always non-negative (additive) or at least $1$
  (multiplicative).

`add_capacity`, `mult_capacity`
: The maximum leakage over a class of priors and/or gain functions.

`strategy(pi, C)`
: An optimal strategy of the adversary: an array with the best action for each observation $y$.

Measures parametrized by a gain or loss function (for instance `g_vuln`) take it as an extra first argument.

```{toctree}
:maxdepth: 1

bayes_vuln
bayes_risk
pred_vuln
pred_risk
g_vuln
l_risk
shannon
guessing
d_privacy
```
