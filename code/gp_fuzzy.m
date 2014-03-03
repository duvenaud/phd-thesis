function gp_fuzzy
%
% A demo showing what happens if you condition on GPs.
%
% David Duvenaud
% March 2014

addpath(genpath( 'utils' ));



% How to save figure.
introfigsdir = '../figures/intro';
savefigs = true;
dpi = 600;

% How figure looks.
num_quantiles = 10;
Nstar = 200;
ylimits = [-1 1.1];
show_samples = true;
sample_ix_subset = false;
num_rand_samples = 10;

% Specify GP model.
sigma = .02;
length_scale = 20;
output_variance = 0.14;

% Make up some data
x = [ 20 80 35 ]';
y = [ -10 55 40 ]' ./100;


max_N = length(x);

% Condition on one datapoint at a time.
for N = 0:max_N

    % Fill in gram matrix
    K = NaN(N,N);
    for j = 1:N
        for k = 1:N
            K(j,k) = kernel( x(j), x(k), length_scale, output_variance );
        end
    end

    % Compute inverse covariance
    Kn = K + sigma^2 .* diag(ones(N, 1));

    % Compute covariance with test points.
    xstar = linspace(-10,110, Nstar);
    Kstar = NaN(Nstar, N);
    kfull = NaN(Nstar, Nstar);
    for j = 1:Nstar
        for k = 1:N
            Kstar(j,k) = kernel( xstar(j), x(k), length_scale, output_variance);
        end
        kstarstar = kernel( xstar(j), xstar(j), length_scale, output_variance);
        for k = 1:Nstar
            kfull(j,k) = kernel( xstar(j), xstar(k), length_scale, output_variance);
        end
    end

    % Compute posterior mean and variance.
    f = (Kstar / Kn) * y(1:N);
    variance = kstarstar - diag((Kstar / Kn) * Kstar');

    full_variance = kfull - (Kstar / Kn) * Kstar';

    % Plot posterior mean and variance.
    figure(N+1); clf;
    quantiles = linspace(0,0.5,num_quantiles+1);%0:0.05:0.5;
    quantiles = quantiles(2:end);

    for s = quantiles
        edges = [f+norminv(s, 0, 1).*sqrt(variance); ...
         flipdim(f-norminv(s, 0, 1).*sqrt(variance),1)]; 
        hc1 = fill([xstar'; flipdim(xstar',1)], edges, color_spectrum(2*s), 'EdgeColor', 'none'); hold on;
    end    






    if show_samples
        % Plot posterior samples
        
        
        seed=0;   % fixing the seed of the random generators
        randn('state',seed);
        rand('state',seed);

        for n_sample = 1:num_rand_samples
            if sample_ix_subset
                cur_ixs = sort(randperm(Nstar, 2));
                cur_range = cur_ixs(1):cur_ixs(2);
            else
                cur_range = 1:Nstar;
            end
            sample = mvnrnd( f, full_variance, 1);
            hs = plot( xstar(cur_range), sample(cur_range), '-', 'Color', colorbrew(n_sample), 'Linewidth', 1); hold on;
            %
        end
    end
    
    h2 = plot( x(1:N), y(1:N), 'kx', 'Linewidth', 1.5, 'Markersize', 10); hold on;

    
    ylim( ylimits);
    xlim( [xstar(1), xstar(end)]);

    set(gca,'Layer','top')   % Show the axes again

    %legend_handle = legend( [hc1 h2 ], { 'GP Posterior Uncertainty', 'Data'}, 'Location', 'SouthEast');
    %set_thesis_fonts( gca, legend_handle );
    set( gca, 'XTick', [] );
    set( gca, 'yTick', [] );
    set( gca, 'XTickLabel', '' );
    set( gca, 'yTickLabel', '' );
    %xlabel( '$x$' );
    %ylabel( '$f(x)$\qquad' );
    set(get(gca,'XLabel'),'Rotation',0,'Interpreter','latex', 'Fontsize', 16);
    set(get(gca,'YLabel'),'Rotation',0,'Interpreter','latex', 'Fontsize', 16);
    set(gcf, 'color', 'white');

    % Add axes, legend, make the plot look nice, and save it.
    tightfig();
    set_fig_units_cm(10,6);

    if savefigs
        filename = sprintf('%s/%s-%d', introfigsdir, 'fuzzy', N );
        %myaa('publish');
        %savepng(gcf, filename);
        save2pdf(filename, gcf, dpi, true );
    end
end
end

function col = color_spectrum(p)
    no_col = [1 1 1];
    full_col = [ 1 0 0 ];
    col = (1 - p)*no_col + p*full_col;
end

function d = kernel(x, y, length_scale, output_variance)
    d = output_variance * exp( - 0.5 * ( (( x - y ) ./ length_scale) .^ 2 )  ); ...
      %+ 3.*output_variance * exp( - 0.5 * ( ( ( x - y ) ./ (length_scale * 10) ).^ 2 ) ) ...
      %+ 0.01 * output_variance * exp( - 0.5 * ( ( ( x - y ) ./ (length_scale ./ 10) ).^ 2 ) );
end
