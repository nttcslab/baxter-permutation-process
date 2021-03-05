# Bayesian Inference for Baxter Permutation Process

![animationMCMCepinions](https://user-images.githubusercontent.com/73105349/96673883-8318c380-13a2-11eb-981f-2f82498ac125.gif)

This is a demo code for Bayesian nonparametric relational data analysis based on **Baxter Permutation Process** ([NeurIPS, 2020](https://nips.cc/)) implemented with Matlab and python.

The key features are listed as follows:

- **Clustering based on rectangular partitioning**: For an input relational matrix, it can discover disjoint rectangle blocks and suitable permutations of rows and columns. 
- **Infinite model complexity**: There is no need to fix the suitable number of rectangle clusters in advance, which is a fundamental principle of Bayesian nonparametric machine learning. 
- **Arbitrary rectangular partitioning**: It can potentially obtain a posterior distribution on arbitrary rectangular partitioning with any numbers of rectangle blocks.  

Please carefully read the licence file before installing and utilizing our software. The following paper must be cited when publishing artiches that adopt or improve out software:

Masahiro Nakano, Akisato Kimura, Takeshi Yamada, Naonori Ueda, "Baxter permutation process," Advances in Neural Information Processing Systems 33 (NeurIPS2020).

For the Matlab implementation, you will need a basic installation without any other Toolboxes. 

For the python implementation, please check `requirements.txt` for the dependency.

## In a nutshell

### Matlab

1. `cd matlab/baxter-permutation-process`
2. `run` 

### python

1. `cd python`
2. `python demo.py`

Then you can see a Markov chain Monte Carlo (MCMC) evolution with the following two figures:
- Rectangular partitioning of a sample matrix (`irmdata\epinions.mat` ).
- Perplexity evolution.

![demoscreen](https://user-images.githubusercontent.com/73105349/96674012-c70bc880-13a2-11eb-921b-3d6b5e78d6cf.gif)

Please check `matlab` or `python` directory for the detailed implementations.

## Reference

Masahiro Nakano, Akisato Kimura, Takeshi Yamada, and Naonori Ueda, 'Baxter Permutation Process,' *Advances in Neural Information Processing Systems* 33 (NeurIPS 2020). 
