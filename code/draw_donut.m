function draw_donut(seed)
%
% Show latent surfaces warped into observed spaces using composite kernels.
%
% David Duvenaud
% March 2014

addpath('utils');

% Fix the seed of the random generators.
if nargin < 1; seed = 4; end
randn('state',seed); rand('state',seed);

topologyfigsdir = '../figures/topology';


output_dim = 3;

exact_grid_resolution = 20;
aux_grid_resulution = 100;
xrange = linspace( 0, 1, exact_grid_resolution);

[x1, x2] = ndgrid( xrange, xrange );   % Make a grid
exactpoints = [ x1(:), x2(:) ];

N = size(exactpoints, 1);

xrange_aux = linspace( 0, 1, aux_grid_resulution );
[x1, x2] = ndgrid( xrange_aux, xrange_aux ); %randn(n, D);   % Deepest latent layer.
x0aux = [ x1(:), x2(:) ];

scale = 8.*pi;
circle_colors = coord_to_color([sin(x0aux(:,1).*scale), cos(x0aux(:,2).*scale)]);
figure(10);
imshow(reshape(circle_colors, [aux_grid_resulution, aux_grid_resulution, 3]))

% Specify the mean and covariance of the warping.
% This mean and variance will be repeated for each dimension.

kernel = @per_kernel_2d;
%kernel = @se_kernel;

mu = zeros(output_dim, N);
sigma = kernel(exactpoints', exactpoints') + eye(N) * 1e-4;
k_x_xaux = kernel(exactpoints', x0aux');


Y = mvnrnd( mu, sigma)';   % Now sample the observed space.
% Work out warping distribution, conditional on the already sampled points.
mucond = k_x_xaux' / sigma * Y;

xaux = mucond;


% Plot the observed space.
figure(11);   % Angry Ivan
tri = delaunay(x0aux(:,1),x0aux(:,2));
trisurf(tri,xaux(:,1),xaux(:,2),xaux(:,3), sum(circle_colors,2), 'EdgeColor', 'none', 'facealpha', 0.5 );

axis off
set(gcf, 'color', 'white');

set_fig_units_cm(10,10);
set(gca, 'projection', 'perspective');
view([-97,4]);
zoom(1.5);
%pos = get(gca, 'CameraPosition')
set(gca, 'CameraPosition', [0.25    0.4    0.5] );
filename = sprintf('%s/%s', topologyfigsdir, 'torus' );
%save2pdf(filename, gcf, 200, true);
myaa('publish');
savepng(gcf, filename);
end


function sigma = se_kernel(x, y)
    if nargin == 0
        sigma = 'Normal SE covariance.'; return;
    end

    sigma = exp( -0.5.*sq_dist(x, y));
end

function sigma = per_kernel(x, y)
    if nargin == 0
        sigma = 'Normal SE covariance.'; return;
    end

    sigma = exp( -0.5.*sin(sq_dist(x, y)).^2);
end

function sigma = per_kernel_2d(x, y)
    if nargin == 0
        sigma = 'Two-dimensional toroidal covariance.'; return;
    end

    period = 1;
    lengthscale = 1000;
    sigma = exp( -2.*sin(pi*(sqrt(sq_dist(x(1,:), y(1,:))))./period).^2 ./ lengthscale)  ...
          .*exp( -2.*sin(pi*(sqrt(sq_dist(x(2,:), y(2,:))))./period).^2 ./ lengthscale);
end
