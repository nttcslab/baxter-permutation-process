function [trueBlockLocation, RowLocation, ColumnLocation, object] = test_MCMC_BPP(X, opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2      % Assign the options to their default values
    opt = [];
end
if ~isfield(opt,'n'),   opt.n = 1;   end;
if ~isfield(opt,'alpha'),   opt.alpha = 0.1;    end;
if ~isfield(opt,'maxiter'), opt.maxiter = 1001; end;
if ~isfield(opt,'ENC'), opt.ENC = size(X,1)*size(X,2)/20000; end;
if ~isfield(opt,'missingRatio'), opt.missingRatio = 0.1; end;
n = opt.n;
alpha = opt.alpha;
maxiter = opt.maxiter;
ENC = opt.ENC;
missingRatio = opt.missingRatio;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% test data %%%%%%%%%
originalX = X;
missingIndex=rand(size(originalX));
X(missingIndex<missingRatio)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Initialization %%%%%%%%%%
BetaRandom = betarnd(ones(n,1), ones(n,1));
BP = [1];
UniformRandom = rand(1);
trueBlockLocation = zeros(1, 4);
trueBlockLocation(1,:) = [0 1 0 1];
for ii=1:n
    if rand(1)>0
        [BP, UniformRandom] = func_EvolvBP(BP, UniformRandom);
    else
        if size(BP,1)>1
            [BP, UniformRandom] = func_DeleteMax(BP, UniformRandom);
        end
    end
    LabelMatrix = func_BP2FP(BP);
    trueBlockLocation = func_EvolvRP(BP, LabelMatrix, BetaRandom(ii), trueBlockLocation);
end
%%%%%%%%% Initialization %%%%%%%%%%

n=size(BP,1);
RowLocation = rand(size(X,1),1);
ColumnLocation = rand(size(X,2),1);

object = [];
for iter=1:maxiter

    [RowLocation, ColumnLocation]...
        = MH_RowLocation(X, trueBlockLocation, RowLocation, ColumnLocation, alpha);

    [BP, BetaRandom, trueBlockLocation, UniformRandom]...
        = MH_RP_demo(X, BP, BetaRandom, trueBlockLocation, UniformRandom, RowLocation, ColumnLocation, alpha, ENC);

    [BP, BetaRandom, trueBlockLocation, UniformRandom]...
        = update_RP(BetaRandom, UniformRandom, X, RowLocation, ColumnLocation, alpha);

    currentPerplexity = calc_testPerplexity(originalX, trueBlockLocation,...
        RowLocation, ColumnLocation, alpha, missingIndex, missingRatio);
    if iter == 1
        object = currentPerplexity;
    else
        object = [object;currentPerplexity];
    end

    subplot(1,2,1);plot_ClusteringResult(originalX, trueBlockLocation, RowLocation, ColumnLocation);
    subplot(1,2,2);plot(object);pause(0.01);
end
