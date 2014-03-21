function decomp_birth(savefigs)
%
% A simple demo script to show how to build a GP model in a realistic scenario.
%
% David Duvenaud
% March 2014

if nargin < 1; savefigs = false; end

randn('state',0);
rand('state',0);

show_samples = false;   % Show samples from the posterior.

% How to save figure.
decompfigsdir = '../figures/worked-example/';
fileprefix = [decompfigsdir 'births'];

dataset_name = 'data/quebec/quebec-all.mat';
load(dataset_name);
feature_names = {'time', 'Weekend', 'Floating holiday' , 'Easter'};

% Fix dates
X(:,1) = linspace(min(X(:,1)),max(X(:,1)), length(X));

X = X(end-2000:end,:);
y = y(end-2000:end);

% Normalize the data.
%X = X - repmat(mean(X), size(X,1), 1 );
%X = X ./ repmat(std(X), size(X,1), 1 );

y = y - mean(y);
%y = y / std(y);

[N,D] = size(X);

% Easy part of model set up:
meanfunc = {'meanZero'};
inference = @infExact;
likfunc = @likGauss;
hyp.mean = [];
hyp.lik = ones(1,eval(likfunc())).*log(0.1);   


longterm_ell = 1.5;  longterm_sf = 14;
longterm_cov = { 'covMask', { 1, 'covSEiso'}};
longterm_hyp = [log(longterm_ell); log(longterm_sf)];

weekly_ell = 1;  weekly_sf = 55; weekly_period = 1/52;
weekly_se_ell = 2;  weekly_se_sf = 1;
%weekly_period_cov = { 'covMask', { 1, 'covPeriodic'}};
weekly_period_cov = { 'covMask', { 1, {'covProd', {'covPeriodic', 'covSEiso'}}}};
weekly_period_hyp = [log(weekly_ell); log(weekly_period); log(weekly_sf); ...
                     log(weekly_se_ell); log(weekly_se_sf)];

yearly_ell = 0.5;  yearly_sf = 17.5; yearly_period = 1;
yearly_se_ell = 5.7;  yearly_se_sf = 1;
yearly_period_cov = { 'covMask', { 1, {'covProd', {'covPeriodic', 'covSEiso'}}}};
yearly_period_hyp = [log(yearly_ell); log(yearly_period); log(yearly_sf); ...
                     log(yearly_se_ell); log(yearly_se_sf)];

%yearly_short_ell = 1/365;  yearly_short_sf = 1;
%yearly_period_short_cov = { 'covMask', { 1, 'covPeriodic'}};
%yearly_period_short_hyp = [log(yearly_short_ell); log(yearly_period); log(yearly_short_sf)];

%medterm_ell = 0.2;  medterm_sf = 10;
%medterm_cov = { 'covMask', { 1, 'covSEiso'}};
%medterm_hyp = [log(medterm_ell); log(medterm_sf)];

shortterm_ell = 0.04;  shortterm_sf = 10;
shortterm_cov = { 'covMask', { 1, 'covSEiso'}};
shortterm_hyp = [log(shortterm_ell); log(shortterm_sf)];



% Now construct the kernel.
kernel_components = cell(0);
kernel_hypers = cell(0);
kernel_names = cell(0);

kernel_names{end+1} = 'Long-term';
kernel_components{end+1} = longterm_cov;
kernel_hypers{end+1} =  longterm_hyp;

kernel_names{end+1} = 'Weekly';
kernel_components{end+1} = weekly_period_cov;
kernel_hypers{end+1} =  weekly_period_hyp;

kernel_names{end+1} = 'Yearly';
kernel_components{end+1} = yearly_period_cov;
kernel_hypers{end+1} =  yearly_period_hyp;

%kernel_components{end+1} = yearly_period_short_cov;
%kernel_hypers{end+1} =  yearly_period_short_hyp;

%kernel_components{end+1} = medterm_cov;
%kernel_hypers{end+1} =  medterm_hyp;

kernel_names{end+1} = 'Short-term';
kernel_components{end+1} = shortterm_cov;
kernel_hypers{end+1} =  shortterm_hyp;

% concatenate all additive components
complete_cov = { 'covSum', kernel_components };

hyp.cov = unwrap(kernel_hypers);

subset = randperm(N, 1000);
%hyp = minimize(hyp, @gp, -100, inference, meanfunc, complete_cov, likfunc, ...
               %X(subset, :), y(subset));

load('hypers.mat', 'hyp');
           
% Pack up hyperparameters again.
kernel_hypers = rewrap(kernel_hypers, hyp.cov);

long_hypers = kernel_hypers{1};
fprintf('Estimated long-term ell: %f\n', exp(long_hypers(1)));
fprintf('Estimated long-term sf: %f\n', exp(long_hypers(2)));
fprintf('\n');

week_hypers = kernel_hypers{2};
fprintf('Estimated weekly ell: %f\n', exp(week_hypers(1)) * exp(week_hypers(5)));
fprintf('Estimated weekly frequency: %f\n', 1/exp(week_hypers(2)));
fprintf('Estimated weekly sf: %f\n', exp(week_hypers(3)));
fprintf('Estimated weekly se ell: %f\n', exp(week_hypers(4)));
fprintf('\n');

year_hypers = kernel_hypers{3};
fprintf('Estimated yearly ell: %f\n', exp(year_hypers(1)) * exp(year_hypers(5)));
fprintf('Estimated yearly period: %f\n', exp(year_hypers(2)));
fprintf('Estimated yearly sf: %f\n', exp(year_hypers(3)));
fprintf('Estimated yearly se ell: %f\n', exp(year_hypers(4)));
fprintf('\n');

%year_period_hypers = kernel_hypers{4};
%fprintf('Estimated yearly period ell: %f\n', exp(year_period_hypers(1)));
%fprintf('Estimated yearly period: %f\n', exp(year_period_hypers(2)));
%fprintf('Estimated yearly period sf: %f\n', exp(year_period_hypers(3)));
%fprintf('\n');

short_hypers = kernel_hypers{4};
fprintf('Estimated short-term ell: %f\n', exp(short_hypers(1)));
fprintf('Estimated short-term sf: %f\n', exp(short_hypers(2)));
fprintf('\n');

plot_additive_decomp_cov( X, y, kernel_components, kernel_hypers, hyp.lik, ...
                      savefigs, fileprefix, kernel_names )

%plot_additive_decomp_fancy( X, y, kernel_components, kernel_hypers, hyp.lik, ...
%                      show_samples, savefigs, fileprefix )
