% A simple demo to help me understand variational methods.
%
% David Duvenaud
% March 2013

% fix the seed of the random number generators.
%seed=0;  randn('state',seed);  rand('state',seed);

% Specify settings.
function gplvm_demo(seed)

if nargin < 1; seed=0; end  % fixing the seed of the random generators
randn('state',seed);
rand('state',seed);

n = 1000;
nxbins = 20;
nybins = 20;

nstar = 300;

histprop = 0.25;
sidemargin = 0.08;
intermargin = 0.04;

% Generate data.
x_prior = @(x) normpdf( x, 0, 1 );
true_x = randn(n,1);
x_lim = [min(true_x) max(true_x)];

n_sample_points = 200;
sample_points = linspace( x_lim(1), x_lim(2), n_sample_points);

% Generate a random function.
numerical_noise = 1e-5;
se_length_scale = 2.5;
se_outout_var = 2;
se_kernel = @(x,y) se_outout_var*exp( - 0.5 * ( ( x - y ) .^ 2 ) ./ se_length_scale^2 );
K = bsxfun(se_kernel, sample_points', sample_points ) + eye(n_sample_points).*numerical_noise; % Evaluate prior.
%Kstar = bsxfun(se_kernel, xstar', xrange ) + eye(n_xstar).*numerical_noise; % Evaluate prior.
y = mvnrnd( zeros(size(sample_points)), K, 1)';
weights = K \ y;  % Now compute kernel function weights.
true_f = @(x)(bsxfun(se_kernel, x', sample_points' )' * weights); % Construct posterior function.


% Condition on samples.
%f = (Kstar / Kn) * y(1:N);

%true_f = @(x) 4 - x.^2;
y = true_f(true_x);

%x_lim = [-2 2];
%y_lim = [0 4];

y_lim = [min(y) max(y)];
xrange = linspace(x_lim(1), x_lim(2), nstar);


figure(100); clf;

h=axes('position',[histprop+sidemargin+intermargin, histprop+sidemargin+intermargin, 1-histprop-sidemargin-2*intermargin, 1-histprop-1*sidemargin-2*intermargin]);
% [left bottom width height]
% where left and bottom define the distance from the lower-left corner of the
% container.

% Plot data
plot(true_x, y_lim(1).*ones(size(true_x)), 'k*' ); hold on;
plot(x_lim(1).*ones(size(y)), y, 'k*' );

gray = [0.9 0.9 0.9];
% Plot connecting lines.
for i = 1:n;
    line([true_x(i), true_x(i)], [y_lim(1), true_f(true_x(i))], 'Color', gray );
    line([x_lim(1), true_x(i)], [true_f(true_x(i)), true_f(true_x(i))], 'Color', gray );
end

% Plot true warping function.
plot(xrange, true_f(xrange'), '-', 'Linewidth', 2); hold on;
%plot(true_x, true_f(true_x), '-', 'Linewidth', 2); hold on;

xlim(x_lim);
ylim(y_lim);
set( gca, 'yTick', [] );
set( gca, 'yTickLabel', '' );
set( gca, 'XTick', [] );
set( gca, 'XTickLabel', '' );

% Plot x dist 
h=axes('position',[histprop+sidemargin+intermargin, sidemargin, 1-histprop-sidemargin-2*intermargin, histprop]);
[xcounts, xbins] = hist(true_x, nxbins);
xprob = xcounts ./ sum(xcounts);
bar(xbins, xprob); hold on;
xb = get(gca,'child');
set(xb,'FaceColor', colorbrew(2));

plot(xrange, x_prior(xrange)./max(x_prior(xrange))*max(xprob), 'r--');

xlim(x_lim);
set( gca, 'yTick', [] );
set( gca, 'yTickLabel', '' );
xlabel('Latent x coordinate');



% Plot y dist
h=axes('position',[sidemargin histprop+sidemargin+intermargin, histprop, 1-histprop-sidemargin-2*intermargin]);
[ycounts, ybins] = hist(y, nybins);
yprob = ycounts ./ sum(ycounts);
h = barh(ybins, yprob); hold on;
set(get(h,'Parent'),'xdir','r')
yb = get(gca,'child');
set(yb,'FaceColor', colorbrew(3));
ylim(y_lim);
set( gca, 'XTick', [] );
set( gca, 'XTickLabel', '' );
ylabel('observed y coordinate');

% Make plot prettier.
set(gcf, 'color', 'white');
set(gca, 'TickDir', 'out');

set_fig_units_cm( 14, 8 );
save2pdf([ 'figures/gplvm_1d_draw_', num2str(seed), '.pdf'], gcf, 600, true);
