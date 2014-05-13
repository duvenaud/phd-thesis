Automatic Model Construction with Gaussian Processes
=====================

<img src="https://raw.githubusercontent.com/duvenaud/phd-thesis/master/figures/topology/mobius.png" width="200"> <img src="https://raw.githubusercontent.com/duvenaud/phd-thesis/master/figures/additive/3d-kernel/3d_add_kernel_321.png" width="200"> <img src="https://raw.githubusercontent.com/duvenaud/phd-thesis/master/figures/deep-limits/map_connected/latent_coord_map_layer_39.png" width="200">



I'm about 10 days from submitting my thesis.  In the spirit of radical openness, I've made the entire repo public.  Think I'm missing something?  Want me to cite you?  <a href="mailto: dkd23@cam.ac.uk">Let me know!</a>  Any feedback would be much appreciated.

Individual chapters:

1. [Introduction](intro.pdf)
2. [Expressing Structure with Kernels](kernels.pdf)
3. [Automatic Model Building](grammar.pdf)
4. [Automatic Model Description](description.pdf)
5. [Deep Gaussian Processes](deeplimits.pdf)
6. [Additive Gaussian Processes](additive.pdf)
7. [Warped Mixture Models](warped.pdf)

Or get the whole thing in [one big PDF](thesis.pdf)



Abstract:
----------

This thesis develops a broad class of models, useful for learning and forecasting in such domains as time series, geological formations, and physical dynamics. These models are based on Gaussian processes, which can express a wide variety of statistical structure, depending on the choice of kernel. Historically, the type of kernel has been chosen by hand by experts. In this thesis, we show how to automate this task, creating an artificial statistician capable of systematically exploring a large space of models.

The introductory chapters show many types of structure, such as periodicity, changepoints, additivity, and symmetries can be encoded by kernels, and that these structure can be combined by combining kernels. For example, compound kernels can produce priors over topological manifolds such as cylinders, toruses, and M\"{o}bius strips, as well as their higher-dimensional analogues.

The main contribution of this thesis is to show how this open-ended space of models can be explored automatically. To do so, we define a simple grammar over kernels, a search criterion (marginal likelihood), and a breadth-first search procedure. Combining these, we present a procedure which takes a dataset, and outputs a detailed report with graphs and automatically-generated text illustrating the qualitatively different, and potentially novel, types of structure discovered in that dataset. This system automates parts of the model-building and analysis currently performed by expert statisticians.

This thesis also explores several extensions to Gaussian process models. First, building on earlier work relating Gaussian processes and neural nets, we explore the natural extensions of these models to *deep kernels* and *deep Gaussian processes*. Second, we examine the model class consisting of the sum of functions of all possible combinations of input variables. We show a close connection between this model class and the recently-developed regularization method of *dropout*. Third, we combine Gaussian processes with the Dirichlet process to produce the *warped mixture model* -- an unsupervised clustering model with nonparametric cluster shapes, and a corresponding latent space in which each cluster has an interpretable parametric form.


Source Code for Experiments
------------------

The experiments for each chapter all live in different github repos:

* Expressing Structure with Kernels: [github.com/duvenaud/phd-thesis/code/](www.github.com/duvenaud/phd-thesis/code/)
* Automatic Model Building: [github.com/jamesrobertlloyd/gp-structure-search/](www.github.com/jamesrobertlloyd/gp-structure-search/)
* Automatic Model Description: [github.com/jamesrobertlloyd/gpss-research/](www.github.com/jamesrobertlloyd/gpss-research/)
* Deep Gaussian Processes: [github.com/duvenaud/deep-limits/](www.github.com/duvenaud/deep-limits/)
* Additive Gaussian Processes: [github.com/duvenaud/additive-gps/](www.github.com/duvenaud/additive-gps/)
* Warped Mixture Models: [github.com/duvenaud/warped-mixtures/](www.github.com/duvenaud/warped-mixtures/)

