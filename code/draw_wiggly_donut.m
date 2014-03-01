function deep_gp_sample
%
% Plot samples and sampled distributions from a deep GP model.
%
% David Duvenaud
% Dec 2102

addpath('utils');

% Fix the seed of the random generators.
seed=2;
randn('state',seed);
rand('state',seed);

n = 300;
n_aux = 1000;
D = 3;
Q = 2;
layers = 1000;

scale = 1.5;
% Deepest latent layer.
%x = randn(n, Q) ./ scale;
%xaux = randn(n_aux, Q) ./ scale;

xrange = linspace( 0, 1, 20);


% Make a grid
[xa, xb] = ndgrid( xrange, xrange ); %randn(n, D);   % Deepest latent layer.
x0 = [ xa(:), xb(:) ];

n = size(x0, 1);

xrange2 = linspace( 0, 1, 50 );
[xa, xb] = ndgrid( xrange2, xrange2 ); %randn(n, D);   % Deepest latent layer.
x0aux = [ xa(:), xb(:) ];

%x0aux = rand(1000, Q) .* 2 - 1;

circle_colors = coord_to_color(x0aux(:,1:2));

%    plot_little_circles(x0aux(:,1), x0aux(:,2), ...
%         0.01, circle_colors, 1 ); hold on; 

%plot_density( x, xaux, 1, circle_colors, 'Original latent space' );

layer = 1;

kernel = @(x,y)( per_kernel_1d(x, y) + per_kernel_2d(x,y));

%for layer = 2:layers
% Specify the mean and covariance of the warping.
% This mean and variance will be repeated for each dimension.
mu = zeros(D, n);
sigma = kernel(x0', x0') + eye(n) * 1e-4;
k_x_xaux = kernel(x0', x0aux');

% Now sample next deepest layer.
for d = 1:D
    x(:,d) = mvnrnd( mu(d, :), sigma);
end

% Work out warping distribution,
% conditional on the already sampled points.
mucond = k_x_xaux' / sigma * x;
%sigmacond = k_x_xaux' / sigma * k_x_xaux;

% Now sample some more points, conditioned on those.
%xaux(:,1) = mvnrnd( mucond(:,1), sigmacond);
%xaux(:,2) = mvnrnd( mucond(:,2), sigmacond);
for d = 1:D
    xaux(:,d) = mucond(:,d);
end

% Plot the warped points and auxiliary points.
%    plot_density( x, xaux, layer, circle_colors );

tri = delaunay(x0aux(:,1),x0aux(:,2));
trimesh(tri,xaux(:,1),xaux(:,2),xaux(:,3))


end


function sigma = se_kernel(x, y)
    if nargin == 0
        sigma = 'Normal SE covariance.'; return;
    end

    sigma = exp( -0.5.*sq_dist(x, y));
end

function sigma = per_kernel_1d(x, y)
    if nargin == 0
        sigma = '1d periodic covariance.'; return;
    end

    period = 0.1;
    lengthscale = 10;
    sigma = 0.0005.*exp( -2.*sin(pi*(sqrt(sq_dist(x(2,:), y(2,:))))./period).^2 ./ lengthscale);
end

function sigma = per_kernel_2d(x, y)
    if nargin == 0
        sigma = 'Toroidal covariance.'; return;
    end

    period = 1;
    lengthscale = 10;
    sigma = exp( -2.*sin(pi*(sqrt(sq_dist(x(1,:), y(1,:))))./period).^2 ./ lengthscale)  ...
          .*exp( -2.*sin(pi*(sqrt(sq_dist(x(2,:), y(2,:))))./period).^2 ./ lengthscale);
end


function plot_density( x, xaux, layer, circle_colors, fig_title )

    if nargin < 5
        fig_title = sprintf('Layer %d', layer);
    end
    
    %figure(layer); clf;
    
    circle_size = .005;
    circle_alpha = 0.16;
    circle_color = colorbrew(2);
    
    %plot3(xaux(:,1), xaux(:,2), xaux(:,3), 'g.');
    %plot_little_circles_3d(x(:,1), x(:,2), x(:,3), ...
    %     circle_size, [0 0 0] ); hold on;    
    %good_xlim = xlim;
    %good_ylim = ylim;
    %good_zlim = zlim;
    
    %plot_little_circles_3d(xaux(:,1), xaux(:,2), xaux(:,3), ...
    %     circle_size, circle_colors ); hold on;    
    hold on;
%    set(h_dots, 'visible', 'off');
    %h_dots = plot(x(:,1), x(:,2), 'k.');
    
    
    %alpha(h_dots, 1);
    
    %axis tight
    axis off
    set( gca, 'XTick', [] );
    set( gca, 'yTick', [] );
    set( gca, 'XTickLabel', '' );
    set( gca, 'yTickLabel', '' );
    set(gcf, 'color', 'white');
    set(gca, 'YGrid', 'off');
    
    %border_scale = 10;
    %xlim( good_xlim );
    %ylim( good_ylim );
    %zlim( good_zlim );
    
    %set_fig_units_cm( 14, 14 );
    %axis tight
    set(gca, 'Units', 'normalized', 'Position', [0 0 1 1])
    %title(fig_title, 'Interpreter', 'Latex', 'FontSize', 18);
    %savepng(gcf, sprintf('figures/deep_gp_sample_layer_%d', layer));
end
