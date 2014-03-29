function additive_extrapolation(savefigs)
% Generate some plots showing additive GP regression on synthetic datasets
% with 1st-order interactions and censored data.
%
% David Duvenaud
% March 2014
% ===============

addpath(genpath( 'utils' ));
addpath(genpath( 'gpml' ));

if nargin < 1; savefigs = false; end

% fixing the seed of the random generators
seed=6;
randn('state',seed);
rand('state',seed);

% Rendering setup
n_1d = 100;      % Fineness of grid.


% Regression problem setup.
n = 100;            % number of observations
s = 2;              % number of relevant variables
noise_std = .02;		% standard deviation of noise
X = rand(n,2) * 4 - 2;
noiseless_Y = truefunc(X);
Y =  noiseless_Y + randn(n,1) * noise_std;   % add some noise with known standard deviation

xlims = [-2.3, 6];
ylims = [-2.3, 2];

% Censor data.
censor_threshold = -1.5;
trainset = or(X(:,1) < censor_threshold, X(:,2) < censor_threshold);
X = X(trainset,:);
y = Y(trainset);
noiseless_Y = noiseless_Y(trainset);

[N,D] = size(X);


% Set up model.
likfunc = 'likGauss'; sn = 0.1; hyp.lik = log(sn);
inference = @infExact;
meanfunc = {'meanConst'}; hyp.mean = 0;

% Train additive model
R = 2;
covfunc = { 'covADD',{1:R,'covSEiso'} };  % Construct an additive kernel
hyp.cov = [ log(ones(1,2*D)), log(ones(1,R))];    % Set hyperparameters.
hyp_add = minimize(hyp, @gp, -20, inference, meanfunc, covfunc, likfunc, X, y);

% Train ARD model
covfunc = { 'covSEard' };
hyp.cov = [ log(ones(1,D)), log(1)];    % Set hyperparameters.
hyp = minimize(hyp, @gp, -20, inference, meanfunc, covfunc, likfunc, X, y);


% generate a grid on which to render to predictive surface.
range = linspace(xlims(1), xlims(2), n_1d);
[a,b] = meshgrid(range, range);
xstar = [ a(:), b(:) ];

% Make ARD predictions.
predictions_ard = gp(hyp, inference, meanfunc, covfunc, likfunc, X, y, xstar);

% Make additive predictions.
covfunc = { 'covADD',{1:R,'covSEiso'} };  % Construct an additive kernel.
predictions_add = gp(hyp_add, inference, meanfunc, covfunc, likfunc, X, y, xstar);

figure(1); clf;
nice_plot_surface(a,b,X, predictions_ard,length( range), xlims, ylims);
hold on; show_xs(X, repmat(ylims(1), size(X,1), 1), noiseless_Y - 0.1);
if savefigs
    nice_figure_save('1st_order_censored_ard')
end

figure(2); clf;
nice_plot_surface(a,b,X, predictions_add, length( range), xlims, ylims);
hold on; show_xs(X, repmat(ylims(1), size(X,1), 1), noiseless_Y - 0.1);
if savefigs
    nice_figure_save('1st_order_censored_add')
end

figure(3); clf;
nice_plot_surface(a,b,X, truefunc(xstar), length( range), xlims, ylims);
hold on;
if savefigs
    nice_figure_save('1st_order_censored_truth')
end


end

function nice_plot_surface(a,b,X,Y,l, xlims, ylims)
    
    surf(a,b,reshape(Y, l, l), 'EdgeColor','none','LineStyle','none','FaceLighting','phong'); 
    xlim(xlims); ylim(xlims); zlim(ylims);
    
    set(gcf, 'color', 'white');
    set(get(gca,'XLabel'),'Rotation',0,'Interpreter','tex', 'Fontsize', 12, 'FontName','Times New Roman');
    set(get(gca,'YLabel'),'Rotation',0,'Interpreter','tex', 'Fontsize', 12, 'FontName','Times New Roman');
    set(get(gca,'ZLabel'),'Rotation',0,'Interpreter','tex', 'Fontsize', 12, 'FontName','Times New Roman');
    set( gca, 'xTickLabel', '' );
    set( gca, 'yTickLabel', '' );    
    set( gca, 'zTickLabel', '' );    
    xlabel('x_1');
    ylabel('x_2');
    zlabel('f(x)  ');
  

    view([-30,42]);
    tightfig
    set_fig_units_cm(6,5);
end

function show_xs(X, lowy, uppery)
    line_width = 2;
    markersize = 14;
    plot3(X(:,1), X(:,2), lowy, 'b.', 'Linewidth', line_width, ...
        'Markersize', markersize);   % show the data
    for i = 1:size(X,1);
        line([X(i,1), X(i,1)], [X(i,2), X(i,2)], [lowy(i), uppery(i)], ...
            'LineStyle', '--', 'Color', [0.5 0.5 0.5]);
    end
end

function y = truefunc(x)
    y = sin(x(:,1) .* 2 - 1.1) + sin(x(:,2) .* 2  - 1.1);
end


function nice_figure_save(filename) 
    myaa('mypublish', ['../figures/additive/', filename]);
    %savepng(gcf, ['../figures/additive/', filename]);
end
