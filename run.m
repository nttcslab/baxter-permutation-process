addpath('BPPCodes');

load('irmdata/epinions.mat');
X = X + 1;
[trueBlockLocation, RowLocation, ColumnLocation, object] = test_MCMC_BPP(X);
