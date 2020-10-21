# Bayesian Inference for Baxter Permutation Process

![animationMCMCepinions](https://user-images.githubusercontent.com/73105349/96673883-8318c380-13a2-11eb-981f-2f82498ac125.gif)

This is a MATLAB code for Bayesian nonparametric relational data analysis based on **Baxter Permutation Process** ([NeurIPS, 2020](https://nips.cc/)). The key features are listed as follows:

- **Clustering based on rectangular partitioning**: For an input relational matrix, it can discover disjoint rectangle blocks and suitable permutations of rows and columns. 
- **Infinite model complexity**: There is no need to fix the suitable number of rectangle clusters in advance, which is a fundamental principle of Bayesian nonparametric machine learning. 
- **Arbitrary rectangular partitioning**: It can potentially obtain a posterior distribution on arbitrary rectangular partitioning with any numbers of rectangle blocks.  

You will need a basic MATLAB installation without any other Toolboxes. 

## Very Quick Guide

1. Unzip to any directory.
2. Go to `\BPPinference` directory.
3. Typing `run`

Then you can see a Markov chain Monte Carlo (MCMC) evolution with the following two figures:
- Rectangular partitioning of a sample matrix.
- Perplexity evolution.

![demoscreen](https://user-images.githubusercontent.com/73105349/96674012-c70bc880-13a2-11eb-921b-3d6b5e78d6cf.gif)

## Inference (Model Fitting) Usage

Given an input relational matrix, the Baxter permutation process can be fitted to it by a MCMC inference algorithm as follows:

`[trueBlockLocation, RowLocation, ColumnLocation, object] = test_MCMC_BPP(X);`

or

`[trueBlockLocation, RowLocation, ColumnLocation, object] = test_MCMC_BPP(X,opt);`

- X: An M by N matrix. Each element must be natural numbers.
- opt.maxiter: Maximum number of iterations. 
- opt.missingRatio: Ratio of test/(training+test). 

Note that the input data must be categorical X(i,j)={1,2,3,...}.  For example, if the input matrix X consists of binary {0,1}, then we can simply rewrite X=X+1 in advance. 

## Analysis Result Visualization

You can easily obtain the visualization of the resulting rectangular partitioning as follows. Our MCMC algorithm has a set of parameters on each MCMC iteration, consisting of `trueBlockLocation, RowLocation, ColumnLocation, object ` .  For example, `training_sample.mat` includes an analysis result sample. 

1. `addpath('BPPCodes');`
2. `load('training_sample.mat');`
3. `plot_ClusteringResult(X, trueBlockLocation, RowLocation, ColumnLocation); daspect([1 1 1]);`

## Sampling from Baxter Permutation Process

Our Baxter permutation process can simulate a sample of a sequence of Baxter permutations (BPs), from the viewpoint of a generative probabilistic model as a prior model.  For example, we draw a BP sample on {1,2, ..., 10} as follows: 

1. `addpath('BPPCodes');`
2. `BP=1; UniformRandom =rand(1);` 
3. `for ii=1:10; [BP, UniformRandom] = func_EvolvBP(BP, UniformRandom); end`
4. Typing `BP`

Then, you can see a BP sample. Moreover, according to the fact that there is one-to-one correspondence between BPs and floorplan partitionings (FPs), the corresponding FP can be also induced from the BP sample:

1. `FP = func_BP2FP(BP);`
2. `imagesc(FP);` 

Finally, you can see a FP sample. 

## Reference

1. M. Nakano, A. Kimura, T. Yamada, and N. Ueda, 'Baxter Permutation Process,' *Advances in Neural Information Processing Systems* 33 (NIPS 2020). 