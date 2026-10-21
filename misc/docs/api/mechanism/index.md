# Mechanism construction

```{py:module} qif.mechanism
```

`qif.mechanism`

Functions for constructing mechanisms: channels $C$ from secrets to `n_cols` observations that are optimal with
respect to some privacy and utility criteria. Most of them are solved by linear programming, and return an empty
array if no mechanism satisfies the constraints.

Utility is measured by a loss function `loss(x, y)`, giving the loss of reporting the observation $y$ when the
secret is $x$. Its expected value $\sum_{x,y} \pi_x C_{x,y}\, \mathrm{loss}(x,y)$ for the prior $\pi$ is the
*expected loss* of the mechanism (see {py:func}`qif.utility.expected_distance`). Note that this is different than
the loss functions of the adversary (see {doc}`../measure/l_risk`). Some functions also take a `hard_max_loss`:
observations with $\mathrm{loss}(x,y)$ larger than it are forbidden, that is $C_{x,y} = 0$.

```{toctree}
:maxdepth: 1

bayes_vuln
bayes_risk
g_vuln
l_risk
shannon
geo_ind
d_privacy
```
