function plot_symmetry_samples
% Simple demo of drawing from a 2-dimensional GP
%
% David Duvenaud
% March 2014
% -=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-

close all

N_1d = 50;   % Fineness of grid in each dimension.

range = linspace( -5, 5, N_1d);   % Choose a set of x locations.


plot_2d_gp_and_kernel( range, @gnn_kernel );
plot_2d_gp_and_kernel( range, @add_kernel );
plot_2d_gp_and_kernel( range, @fa_kernel );
plot_2d_gp_and_kernel( range, @se_kernel );
plot_2d_gp_and_kernel( range, @se_1d_kernel );
plot_2d_gp_and_kernel( range, @se_1d_kernel_plus_2d );
plot_2d_gp_and_kernel( range, @symm_xy_kernel_naive );
plot_2d_gp_and_kernel( range, @symm_xy_kernel_proj );
plot_2d_gp_and_kernel( range, @spline_kernel );

end




function plot_2d_gp_and_kernel( range, kernel )

    figdir = '../figures/symmetries/';
    titlestring = kernel();  % Call with no arguments to get a description of the kernel.
    filename = [figdir titlestring]
    
    % Set the random seed.
    randn('state', 0);
    rand('twister', 0);  

    [x1, x2] = meshgrid( range);
    x = [x1(:) x2(:)]';
    N = size(x, 2);
    mu = zeros(N, 1);   % Set the mean of the GP to zero everywhere.

    figure;
    
    % Plot kernel.
    %subplot(1,2,1);
    sigma = kernel( x, [ 0; 0 ] );
    surf(x1, x2, reshape(sigma, sqrt(numel(sigma)), sqrt(numel(sigma))));
    
    nice_figure_save([filename, '-kernel'])
    
    %title('kernel');
    % Plot a draw from a GP with that kernel.
    %subplot(1,2,2);
    sigma = kernel( x, x );
    sigma = 0.5*(sigma + sigma');
    sigma = sigma + eye(length(sigma)) .* 1e-6;
    f = mvnrnd( mu, sigma ); % Draw a sample from a multivariate Gaussian.
    surf(x1, x2, reshape(f, sqrt(numel(f)), sqrt(numel(f))));              % Plot the drawn value at each location.
    %title('sample');
    nice_figure_save([filename, '-sample'])
    
    %mtit(titlestring);
    %set(gcf,'units','centimeters')
    %set(gcf,'Position',[1 1 40 15])
    %savepng(gcf, titlestring );
end






function nice_figure_save(filename)

%    axis off
    set(gcf, 'color', 'white');
    set( gca, 'xTickLabel', '' );
    set( gca, 'yTickLabel', '' );    
    set( gca, 'zTickLabel', '' );    
    set(gca, 'TickDir', 'in')

    tightfig
    set_fig_units_cm(12,12);
    
    myaa('publish');
    savepng(gcf, filename);    
    %filename_eps = ['../figures/additive/3d-kernel/', filename, '.eps']
    %filename_pdf = ['../figures/additive/3d-kernel/', filename, '.pdf']
    %print -depsc2 filename_eps
    %eps2pdf( filename_eps, filename_pdf, true);
end


% Stationary covariance functions
% ===================================


function sigma = gnn_kernel(x, y)
    if nargin == 0
        sigma = 'gaussian-nn'; return;
    end
    
    [dx, nx] = size(x);
    [dy, ny] = size(y);
    n = min(nx, ny);
    sigma = exp( -0.15.*sq_dist(x, zeros(dx, n))) ...
         .* exp( -0.15.*sq_dist(y, zeros(dy, n)))' ...  % Note the transpose!
         .* exp( -0.5.*sq_dist(x, y));
end


function sigma = symm_xy_kernel_naive(x, y)
    if nargin == 0
        sigma = 'symmetric-xy-naive'; return;
    end
    
    xs = [x(2, :); x(1, :)];
    sigma = 0.5.*exp( -0.5.*sq_dist(x, y)) + 0.5.*exp( -0.5.*sq_dist(xs, y));
end


function sigma = symm_xy_kernel_proj(x, y)
    if nargin == 0
        sigma = 'symmetric-xy-projection'; return;
    end
    
    projx = [(x(1,:) + x(2,:))./2; abs(x(1,:) - x(2,:))./2];
    projy = [(y(1,:) + y(2,:))./2; abs(y(1,:) - y(2,:))./2];
    sigma = exp( -0.5.*sq_dist(projx, projy));
end

function sigma = se_kernel(x, y)
    if nargin == 0
        sigma = 'squared-exp'; return;
    end

    sigma = exp( -0.5.*sq_dist(x, y));
end

function sigma = se_1d_kernel(x, y)
    if nargin == 0
        sigma = '1d-squared-exp'; return;
    end
    
    A = [ 0.5 -0.7 ];
    sigma = exp( -0.5.*sq_dist(A*x, A*y));
end

function sigma = se_1d_kernel_plus_2d(x, y)
    if nargin == 0
        sigma = '1d-squared-exp-plus-2d'; return;
    end

    A = [ .7 -.9 ];
    sigma = exp( -0.5.*(sq_dist(A*x, A*y)) + 0.1.*exp( -0.5.*sq_dist(x, y)));
end

function sigma = fa_kernel(x, y)
    if nargin == 0
        sigma = 'factor-analysis'; return;
    end
    
    A = [ .7 -.9 ];
    sigma = exp( -0.5.*(sq_dist(A*x, A*y) + 0.05.*sq_dist(x, y)));
end

function sigma = add_kernel(x, y)
    if nargin == 0
        sigma = 'factor-analysis-plus-se'; return;
    end
    
    sigma = 0.4.*exp( -0.5.*( sq_dist(x(1,:), y(1,:)))) + 0.4.*exp( -0.5.*(sq_dist(x(2,:), y(2,:)))) + 0.1.*exp( -0.5.*(sq_dist(x, y)));
end

function sigma = spline_kernel(x, y)
    if nargin == 0
        sigma = 'exp-y-spline'; return;
    end
    
    sigma = exp( -sqrt(sq_dist(x(1,:), y(1,:))/100) -sqrt(sq_dist(x(2,:), y(2,:))/100));
end

