function draw_shapes()
%
% Show latent surfaces warped into observed spaces using composite kernels.
%
% David Duvenaud
% March 2014

addpath('utils');


% Specify coordinates for different shapes
torus.name = 'torus';
torus.kernel = @per_kernel_2d;
torus.seed = 4;
torus.color_scale = 8.*pi;
torus.camera = [0.25 0.4 0.5];
torus.view = [-100,9];

manifold.name = 'manifold';
manifold.kernel = @se_kernel;
manifold.seed = 0;
manifold.color_scale = 8.*pi;
manifold.camera = [-6.2031  -13.0451    1.7971];
manifold.view = [-57,-5];

cylinder.name = 'cylinder';
cylinder.kernel = @(x,y)(se_kernel(x(1,:), y(1,:)).*per_kernel(x(2,:),y(2,:)));
cylinder.seed = 0;
cylinder.color_scale = 8.*pi;
cylinder.camera = [-8.4894    8.3373   12.5418];
cylinder.view = [-134,31];

mobius.name = 'mobius';
mobius.kernel = @mobius_kernel;
mobius.seed = 0;
mobius.color_scale = 8.*pi;
mobius.camera = [-0.0533    0.0835    0.0333];
mobius.view = [75,-9];


draw_shape(cylinder);
draw_shape(manifold);
draw_shape(torus);
draw_shape(mobius);

end



function draw_shape(shape)

% Fix the seed of the random generators.
randn('state', shape.seed);
rand('state', shape.seed);

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

% Specify color pattern for traingles.
circle_colors = coord_to_color([sin(x0aux(:,1).*shape.color_scale), ...
                                cos(x0aux(:,2).*shape.color_scale)]);

mu = zeros(output_dim, N);
sigma = shape.kernel(exactpoints', exactpoints') + eye(N) * 1e-4;
k_x_xaux = shape.kernel(exactpoints', x0aux');


Y = mvnrnd( mu, sigma)';   % Sample the observed space.
% Work out warping distribution, conditional on the already sampled points.
xaux = k_x_xaux' / sigma * Y;

% Plot the observed space.
figure; clf;
tri = delaunay(x0aux(:,1),x0aux(:,2));
trisurf(tri,xaux(:,1),xaux(:,2),xaux(:,3), sum(circle_colors,2), 'EdgeColor', 'none', 'facealpha', 0.5 );

axis off
set(gcf, 'color', 'white');

tightfig
set_fig_units_cm(10,10);
%set(gca, 'projection', 'perspective');
view(shape.view);
set(gca, 'CameraPosition', shape.camera );
xlim([min(xaux(:,1)), max(xaux(:,1))]);
ylim([min(xaux(:,2)), max(xaux(:,2))]);
zlim([min(xaux(:,3)), max(xaux(:,3))]);
%pos = get(gca, 'CameraPosition')

% Render to image file.
filename = sprintf('%s/%s', topologyfigsdir, shape.name );
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
        sigma = 'Periodic covariance.'; return;
    end
    
    period = 1;
    lengthscale = 10;

    sigma = exp( -2.*sin(pi*(sqrt(sq_dist(x, y)))./period).^2 ./ lengthscale);
end

function sigma = per_kernel_2d(x, y, period, lengthscale)
    if nargin == 0; sigma = 'Two-dimensional toroidal covariance.'; return; end
    if nargin < 3; period = 1; end
    if nargin < 4; lengthscale = 1000; end;

    sigma = exp( -2.*sin(pi*(sqrt(sq_dist(x(1,:), y(1,:))))./period).^2 ./ lengthscale)  ...
          .*exp( -2.*sin(pi*(sqrt(sq_dist(x(2,:), y(2,:))))./period).^2 ./ lengthscale);
end


function sigma = mobius_kernel(x, y, period, lengthscale)
    if nargin == 0; sigma = 'Two-dimensional mobius covariance.'; return; end
    if nargin < 3; period = 1; end
    if nargin < 4; lengthscale = 50; end;

    sigma = per_kernel_2d(x, y, period, lengthscale) ...
          + per_kernel_2d([x(2,:); x(1,:)], y, period, lengthscale);
end


function sigma = cone_kernel(x, y)
    if nargin == 0
        sigma = 'Toroidal covariance.'; return;
    end

    nx = size(x, 1);
    ny = size(y, 1);
    
    xmat = repmat(x(:,2), 1, ny);
    ymat = repmat(y(:,2), 1, nx);
    period = 1;
    lengthscale = 10;
    sigma = exp( -2.*sin(pi*(sqrt(sq_dist(x(:,1)', y(:,1)')))./period).^2 ./ lengthscale)  ...
          .*xmat.*ymat';
end
