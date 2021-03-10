# Bayesian Inference for Baxter Permutation Process

![animationMCMCepinions](https://user-images.githubusercontent.com/73105349/96673883-8318c380-13a2-11eb-981f-2f82498ac125.gif)

This is a python code for Bayesian nonparametric relational data analysis based on **Baxter Permutation Process** ([NeurIPS, 2020](https://nips.cc/)). The key features are listed as follows:

- **Clustering based on rectangular partitioning**: For an input relational matrix, it can discover disjoint rectangle blocks and suitable permutations of rows and columns.
- **Infinite model complexity**: There is no need to fix the suitable number of rectangle clusters in advance, which is a fundamental principle of Bayesian nonparametric machine learning.
- **Arbitrary rectangular partitioning**: It can potentially obtain a posterior distribution on arbitrary rectangular partitioning with any numbers of rectangle blocks.  

Please check requirements.txt for the requirements.

## In a nutshell

1. python demo.py

Then you can see a Markov chain Monte Carlo (MCMC) evolution with the following two figures:
- Rectangular partitioning of a sample matrix (`irmdata\epinions.mat` ).
- Perplexity evolution.

![demoscreen](https://user-images.githubusercontent.com/73105349/96674012-c70bc880-13a2-11eb-921b-3d6b5e78d6cf.gif)

## Reference

1. M. Nakano, A. Kimura, T. Yamada, and N. Ueda, 'Baxter Permutation Process,' *Advances in Neural Information Processing Systems* 33 (NeurIPS 2020).
