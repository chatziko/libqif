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

      The solver to use: ``"INTERNAL"`` (a simplex implementation, which also supports ``rat``), ``"GLPK"``,
      ``"GLOP"``, ``"CLP"``, ``"HIGHS"``, ``"GUROBI"`` or ``"CPLEX"`` (the last five via OR-Tools; ``"GUROBI"``
      needs Gurobi to be installed, ``"CPLEX"`` needs OR-Tools to be compiled with CPLEX support). Which ones are
      available depends on how qif was built; the packages on PyPI include OR-Tools (with GLOP, CLP and HiGHS) but
      not GLPK. With ``"AUTO"``, ``rat`` programs use ``"INTERNAL"``; ``double`` ones use ``"CLP"`` (or ``"GLPK"``,
      ``"INTERNAL"`` if OR-Tools is not available), and ``"HIGHS"`` for the ``"INTERIOR"`` method.

   .. attribute:: method
      :value: "AUTO"

      The method to use: ``"SIMPLEX_PRIMAL"``, ``"SIMPLEX_DUAL"``, ``"INTERIOR"`` (interior point), or ``"AUTO"``
      (chosen by the solver). The ``"INTERNAL"`` solver only supports ``"SIMPLEX_PRIMAL"``.

   .. attribute:: presolve
      :value: True

      Whether to use the solver's presolver.

   .. attribute:: msg_level
      :value: "OFF"

      The amount of solver output: ``"OFF"``, ``"ERR"``, ``"ON"`` or ``"ALL"``.
```
