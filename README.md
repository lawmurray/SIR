LibBi package: SIR
==================

Synopsis
--------

    ./run.sh
    
This samples from the posterior distribution using a bridge particle filter.

    octave --path oct/ --eval "plot_and_print"

This plots the results.


Description
-----------

This package includes a stochastic SIR (susceptible/infectious/recovered)
epidemiological compartmental model of the form

\begin{eqnarray}
\end{eqnarray}

It also includes ab observational data set of an epidemic of Russian influenza
at a boys boarding school (Anonymous 1978). As this is a closed system the
observations are considered exact, and the task is to simulate diffusion
bridges between the observed values, and to estimate parameters.

The model and data set were used as a test case in Del Moral & Murray
(2014). The package may be used to reproduce the results in that paper.


References
----------

Anonymous. Influenza in a boarding school. *British Medical Journal*, 1978, 1,
587.

Del Moral, P. & Murray, L. M. Sequential Monte Carlo with Highly Informative
Observations, 2014.
