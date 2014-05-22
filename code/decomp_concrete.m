function decomp_concrete(savefigs)
%
% A simple demo script to show the decomposition of a dataset into individual components.
%
% David Duvenaud
% 2011

if nargin < 1; savefigs = false; end

show_samples = true;

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
hyp.cov = unwrap(kernel_hypers);

hyp = minimize(hyp, @gp, -100, inference, meanfunc, complete_cov, likfunc, X, y);

% Pack up hyperparameters again.
kernel_hypers = rewrap(kernel_hypers, hyp.cov);


plot_additive_decomp( X, y, kernel_components, kernel_hypers, hyp.lik, ...
                      show_samples, savefigs, fileprefix )

plot_cov = false;
if plot_cov
    % Remove coarse and fine
    remove = [6,7];
    kernel_components(remove) = [];
    kernel_hypers(remove) = [];
    feature_names(remove) = [];
    plot_additive_decomp_cov( X, y, kernel_components, kernel_hypers, hyp.lik, ...
                          savefigs, fileprefix, feature_names )
end