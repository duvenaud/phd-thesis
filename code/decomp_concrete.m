%
% A simple demo script to show the decomposition of a dataset into individual components.
%
% David Duvenaud
% 2011


savefigs = true;

% How to save figure.
decompfigsdir = '../figures/decomp/';
fileprefix = [decompfigsdir 'concrete'];

dataset_name = 'data/r_concrete_500.mat';
load(dataset_name);
feature_names = {'Cement','Slag','Fly Ash','Water','Plasticizer','Coarse','Fine','Age'};

X = X(1:100,:);
y = y(1:100);

% Normalize the data.
X = X - repmat(mean(X), size(X,1), 1 );
X = X ./ repmat(std(X), size(X,1), 1 );

y = y - mean(y);
y = y / std(y);

[N,D] = size(X);

% Easy part of model set up:
meanfunc = {'meanZero'};
inference = @infExact;
likfunc = @likGauss;
hyp.mean = [];
hyp.lik = ones(1,eval(likfunc())).*log(0.1);   


ell = 1;  sf = 0.1;

% Now construct the kernel.
kernel_components = cell(1, D);
kernel_hypers = cell(1, D);
for i = 1:D
    kernel_components{i} = { 'covMask', { i, 'covSEiso'}};
    kernel_hypers{i} = [log(ell);log(sf)];
end

% concatenate all additive components
complete_cov = { 'covSum', kernel_components };
hyp.cov = horzcat( kernel_hypers{:} );

hyp = minimize(hyp, @gp, -100, inference, meanfunc, complete_cov, likfunc, X, y);

% Break up hyperparameters again.
for i = 1:D
    kernel_hypers{i} = [hyp.cov(2*D - 1), hyp.cov(2*D)];
end



plot_additive_decomp( X, y, kernel_components, kernel_hypers, hyp.lik, savefigs, fileprefix )
