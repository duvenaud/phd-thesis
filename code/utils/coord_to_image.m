function colors = coord_to_image(x)
% Takes a set of 2d points and turns them into colors
%
% David Duvenaud
% March 2014

x = mod(x,1);

arrow = [...
[ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]; ...
[ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]; ...
[ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]; ...
[ 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 ]; ...
[ 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 ]; ...
[ 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 ]; ...
[ 0 0 0 0 1 1 0 1 1 0 1 1 0 0 0 0 ]; ...
[ 0 0 0 1 1 0 0 1 1 0 0 1 1 0 0 0 ]; ...
[ 0 0 1 1 0 0 0 1 1 0 0 0 1 1 0 0 ]; ...
[ 0 1 1 0 0 0 0 1 1 0 0 0 0 1 1 0 ]; ...
[ 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 ]; ...
[ 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 ]; ...
[ 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 ]; ...
[ 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 ]; ...
[ 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 ]; ...
[ 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 ]; ...
[ 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 ]; ...
[ 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 ]; ...
[ 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 ]; ...
];

% add some padding
arrow = [zeros( size(arrow, 1), 4), arrow, zeros( size(arrow, 1), 4)];
arrow = [zeros( 4, size(arrow, 2)); arrow; zeros( 4, size(arrow, 2))];

[n,d] = size(x);
assert(d == 2);

grid_resolution = 14;
%coords = [[sin(theta) 0]; [cos(theta) 0]];


xrange1 = linspace( 0, 1, size(arrow,1));
xrange2 = linspace( 0, 1, size(arrow,2));
[x1, x2] = ndgrid( xrange1, xrange2 );   % Make a grid
coords = [ x1(:), x2(:) ];
N = size(coords, 1);

corner_colors = repmat(colorbrew(1), N, 1);
corner_colors(logical(arrow(:)), :) = repmat(colorbrew(2), sum(arrow(:)), 1);

%coords = [[xmin ymin]; [xmin ymax]; [xmax ymin]; [xmax ymax]];
weights = se_kernel( x', coords');


% Take each point, and its color will be a linear combination of the 4
% original colors, depending on their location.

weight_sums = sum(weights, 2);
weights = weights ./ repmat( weight_sums, 1, length(coords));

colors(:,1) = weights * corner_colors(:, 1);
colors(:,2) = weights * corner_colors(:, 2);
colors(:,3) = weights * corner_colors(:, 3);
end

function sigma = se_kernel(x, y)
    if nargin == 0
        sigma = 'Normal SE covariance.'; return;
    end

    sigma = exp( -0.5.*sq_dist(x, y).*1000);
end
