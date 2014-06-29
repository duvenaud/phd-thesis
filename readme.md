Automatic Model Construction with Gaussian Processes
=====================

<img src="https://raw.githubusercontent.com/duvenaud/phd-thesis/master/figures/topology/mobius.png" width="200"> <img src="https://raw.githubusercontent.com/duvenaud/phd-thesis/master/figures/additive/3d-kernel/3d_add_kernel_321.png" width="200"> <img src="https://raw.githubusercontent.com/duvenaud/phd-thesis/master/figures/deep-limits/map_connected/latent_coord_map_layer_39.png" width="200">

Defended on June 26th, 2014.

Individual chapters:

1. [Introduction](intro.pdf)
2. [Expressing Structure with Kernels](kernels.pdf)
3. [Automatic Model Building](grammar.pdf)
4. [Automatic Model Description](description.pdf)
5. [Deep Gaussian Processes](deeplimits.pdf)
6. [Additive Gaussian Processes](additive.pdf)
7. [Warped Mixture Models](warped.pdf)
7. [Discussion](discussion.pdf)

Or get the whole thing in [one big PDF](thesis.pdf)



Abstract:
----------

This thesis develops a method for automatically constructing, visualizing and describing a large class of models, useful for forecasting and finding structure in domains such as time series, geological formations, and physical dynamics.
These models, based on Gaussian processes, can capture many types of statistical structure, such as periodicity, changepoints, additivity, and symmetries.
Such structure can be encoded through *kernels*, which have historically been hand-chosen by experts.
We show how to automate this task, creating a system that explores an open-ended space of models and reports the structures discovered.

To automatically construct Gaussian process models, we search over sums and products of kernels, maximizing the approximate marginal likelihood.
We show how any model in this class can be automatically decomposed into qualitatively different parts, and how each component can be visualized and described through text.
We combine these results into a procedure that, given a dataset, automatically constructs a model along with a detailed report containing plots and generated text that illustrate the structure discovered in the data.

The introductory chapters contain a tutorial showing how to express many types of structure through kernels, and how adding and multiplying different kernels combines their properties.
Examples also show how symmetric kernels can produce priors over topological manifolds such as cylinders, toruses, and Mobius strips, as well as their higher-dimensional generalizations.

This thesis also explores several extensions to Gaussian process models.
First, building on existing work that relates Gaussian processes and neural nets, we analyze natural extensions of these models to *deep kernels* and *deep Gaussian processes*.
Second, we examine *additive Gaussian processes*, showing their relation to the regularization method of *dropout*.
Third, we combine Gaussian processes with the Dirichlet process to produce the *warped mixture model*: a Bayesian clustering model having nonparametric cluster shapes, and a corresponding latent space in which each cluster has an interpretable parametric form.


Source Code for Experiments
------------------

The experiments for each chapter all live in different github repos:

* Expressing Structure with Kernels: [github.com/duvenaud/phd-thesis/code/](http://github.com/duvenaud/phd-thesis/tree/master/code/)
* Automatic Model Building: [github.com/jamesrobertlloyd/gp-structure-search/](http://www.github.com/jamesrobertlloyd/gp-structure-search/)
* Automatic Model Description: [github.com/jamesrobertlloyd/gpss-research/](http://www.github.com/jamesrobertlloyd/gpss-research/)
* Deep Gaussian Processes: [github.com/duvenaud/deep-limits/](http://www.github.com/duvenaud/deep-limits/)
* Additive Gaussian Processes: [github.com/duvenaud/additive-gps/](http://www.github.com/duvenaud/additive-gps/)
* Warped Mixture Models: [github.com/duvenaud/warped-mixtures/](http://www.github.com/duvenaud/warped-mixtures/)

