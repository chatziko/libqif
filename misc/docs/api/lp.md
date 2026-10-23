# Linear solver

```{py:module} qif.lp
```

`qif.lp`

Many functions of qif (for instance in {doc}`mechanism/index`, {doc}`refinement` and {doc}`metric/optimize`) solve
linear programs. Their solver can be configured globally, by assigning to the attributes of {py:class}`defaults`,
for instance `qif.lp.defaults.solver = "CLP"`.

## Settings

```{eval-rst}
.. autoclass:: defaults

   .. attribute:: solver
      :value: "AUTO"

      The solver to use: ``"SOPLEX"`` (`SoPlex <https://soplex.zib.de/>`_, which solves ``rat`` programs exactly and
      ``double`` ones in floating point, always available), ``"INTERNAL"`` (a simple simplex implementation, which also
      supports ``rat``), ``"GLOP"``, ``"CLP"``, ``"HIGHS"``, ``"GUROBI"`` or ``"CPLEX"`` (the last five via OR-Tools;
      ``"GUROBI"`` needs Gurobi to be installed, ``"CPLEX"`` needs OR-Tools to be compiled with CPLEX support). The
      packages on PyPI include OR-Tools (with GLOP, CLP and HiGHS). With ``"AUTO"``, ``rat`` programs use
      ``"SOPLEX"``; ``double`` ones use ``"CLP"`` (``"SOPLEX"`` if OR-Tools is not available), and ``"HIGHS"`` for the
      ``"INTERIOR"`` method (which needs OR-Tools). For badly-scaled programs (eg d-privacy with large epsilon),
      ``"SOPLEX"`` is more robust than the OR-Tools solvers.

   .. attribute:: method
      :value: "AUTO"

      The method to use: ``"SIMPLEX_PRIMAL"``, ``"SIMPLEX_DUAL"``, ``"INTERIOR"`` (interior point), or ``"AUTO"``
      (chosen by the solver). The ``"INTERNAL"`` solver only supports ``"SIMPLEX_PRIMAL"``, and ``"SOPLEX"`` does
      not support ``"INTERIOR"``.

   .. attribute:: presolve
      :value: True

      Whether to use the solver's presolver.

   .. attribute:: msg_level
      :value: "OFF"

      The amount of solver output: ``"OFF"``, ``"ERR"``, ``"ON"`` or ``"ALL"``.
```
