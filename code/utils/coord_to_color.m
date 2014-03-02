function colors = coord_to_color(x, scale)
% Takes a set of 2d points and turns them into colors
%
% David Duvenaud
% March 2014


if nargin < 2; scale = 5; end

corner_colors = colorbrew([5 2 3 4 1]);
[n,d] = size(x);
assert(d == 2);

% Rescale points from 0 to 1 along each dimension.

%x(:,1) = x(:,1) - min(x(:,1));
%x(:,2) = x(:,2) - min(x(:,2));

%xmax = max(x(:,1));
%ymax = max(x(:,2));

%x(:, 1) = x(:, 1) ./ xmax;
%x(:, 2) = x(:, 2) ./ ymax;

%xmin = 0.1;
%xmax = 0.9;
%ymin = 0.1;
%ymax = 0.9;

theta = (0:(2*pi/4):2*pi) + pi/4;
theta = theta(1:end-1);
coords = [[sin(theta) 0]; [cos(theta) 0]];

%coords = [[xmin ymin]; [xmin ymax]; [xmax ymin]; [xmax ymax]];
weights = se_kernel( x' .*scale, coords.*scale);


% Take each point, and its color will be a linear combination of the 4
% original colors, depending on their location.

weight_sums = sum(weights, 2);
weights = weights ./ repmat( weight_sums, 1, length(coords));

colors(:,1) = weights * corner_colors(:, 1);
colors(:,2) = weights * corner_colors(:, 2);
colors(:,3) = weights * corner_colors(:, 3);

%{
xcolors = bsxfun( @times, x(:,1), corner_colors(1,:)) ...
        + bsxfun( @times, 1 - x(:,1), corner_colors(2,:));
ycolors = bsxfun( @times, x(:,2), corner_colors(3,:)) ...
        + bsxfun( @times, 1 - x(:,2), corner_colors(4,:));
    
xcolors = (tanh((xcolors - 0.5) .* 2) + 1) ./ 2;    
ycolors = (tanh((ycolors - 0.5) .* 2) + 1) ./ 2;    

colors = (xcolors + ycolors)./2;
%colors = colors .^ 2;
%colors = (tanh((colors - 0.5) .* 5) + 1) ./ 2;
%}
end

function sigma = se_kernel(x, y)
    if nargin == 0
        sigma = 'Normal SE covariance.'; return;
    end

    sigma = exp( -0.5.*sq_dist(x, y));
end