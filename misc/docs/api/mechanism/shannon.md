# Mechanism construction for Shannon entropy

```{py:module} qif.mechanism.shannon
```

`qif.mechanism.shannon`

## Functions

```{eval-rst}
.. autofunction:: max_entropy_given_same_loss

   A channel minimizing the mutual information :math:`I(X;Y)` (that is, maximizing the posterior Shannon entropy) among
   all channels with at most the same expected loss as itself: a point of the rate-distortion curve, computed with
   the Blahut–Arimoto algorithm [CT06]. The trade-off is controlled by scaling ``loss`` (larger values give lower
   expected loss and higher leakage). ``out`` is the initial output distribution; the algorithm iterates until it
   is stable up to ``md``, ``mrd``. Only ``double`` is supported.
```

## References

[CT06]
: T. M. Cover, J. A. Thomas. *Elements of Information Theory*, 2nd edition, Section 10.8. Wiley, 2006.
