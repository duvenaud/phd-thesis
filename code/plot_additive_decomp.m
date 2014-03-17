function plot_additive_decomp( X, y, kernel_components, kernel_params, log_noise, savefigs, fileprefix )
%
% Decomposes an additive GP model into component parts
%
% X,y: The original data
% kernel_components: a cell array of covariance functions, to be added
% together.  Must include noise!
% kernel_params: a cell array of kernel hyperparameters.
%
% David Duvenaud
% March 2014

addpath(genpath( 'utils' ));
addpath(genpath( 'gpml' ));

if nargin < 1; savefigs = false; end

nstar = 200;
show_samples = false;

[N,D] = size(X);
num_components = length(kernel_components);
assert(num_components == length(kernel_params));
assert(D == num_components);

noise_var = exp(2*log_noise);

% First, compute and plot the entire model.
% concatenate all additive components
complete_cov = { 'covSum', kernel_components };
complete_hypers = horzcat( kernel_params{:} );

complete_sigma = feval(complete_cov{:}, complete_hypers, X) + eye(length(y)).*noise_var;

%complete_sigmastar = feval(complete_cov{:}, complete_hypers, X, xrange);
%complete_sigmastarstart = feval(complete_cov{:}, complete_hypers, xrange);

% First, plot the data
%complete_mean = complete_sigmastar' / complete_sigma * y;
%complete_var = complete_sigmastarstart - complete_sigmastar' / complete_sigma * complete_sigmastar;



%figure(1); clf;
%filename = sptrintf( '%s-complete', fileprefix );
%nice_oned_plot( X, y, xstar, complete_mean, complete_var, true, filename, savefigs)



% Next, show the posterior of each component, one at a time.
for d = 1:D
    
    left_extend = 0.1;
    right_extend = 0.1;
    x_left = min(X(:,d)) - (max(X(:,d)) - min(X(:,d)))*left_extend;
    x_right = max(X(:,d)) + (max(X(:,d)) - min(X(:,d)))*right_extend;
    xrange = linspace(x_left, x_right, nstar)';
    xstar = NaN(nstar,D);
    xstar(:, d) = xrange;
    
    cur_cov = kernel_components{d};
    cur_hyp = kernel_params{d};

    % Compute posterior of just this component.
    decomp_sigma = feval(cur_cov{:}, cur_hyp, X);
    decomp_sigma_star = feval(cur_cov{:}, cur_hyp, xstar, X);
    decomp_sigma_starstar = feval(cur_cov{:}, cur_hyp, xstar);
    decomp_mean = decomp_sigma_star / complete_sigma * y;
    decomp_var = decomp_sigma_starstar - decomp_sigma_star / complete_sigma * decomp_sigma_star';
    
    data_mean = decomp_sigma' / complete_sigma * y;


    % Plot posterior mean and variance.
    figure(d); clf;
    filename = sprintf( '%s-component-%d', fileprefix, d );
    nice_oned_plot( X(:,d), y, data_mean, xrange, decomp_mean, decomp_var, show_samples, savefigs, filename)
end



end


function nice_oned_plot( X, y, y_adjusted, xstar, mean, full_variance, show_samples, savefigs, filename)
% Makes a nice plot of a GP posterior.
%
% David Duvenaud
% March 2014


    variance = diag(full_variance);

    % How figure looks.
    num_quantiles = 10;
    num_rand_samples = 10;
    dpi = 600;

    if nargin < 7; savefigs = false; end
    
    
    quantiles = linspace(0,0.5,num_quantiles+1);
    quantiles = quantiles(2:end);

    for s = quantiles
        edges = [mean + norminv(s, 0, 1).*sqrt(variance); ...
         flipdim(mean - norminv(s, 0, 1).*sqrt(variance),1)]; 
        h_gp_post = fill([xstar; flipdim(xstar,1)], edges, color_spectrum(2*s), ...
                   'EdgeColor', 'none'); hold on;
    end    
    
    h_data_orig = plot( X, y, 'ko', 'Linewidth', 1.5, 'Markersize', 10, 'Color', colorbrew(2)); hold on;

    if show_samples
        % Plot posterior samples
        
        
        seed=0;   % fixing the seed of the random generators
        randn('state',seed);
        rand('state',seed);

        for n_sample = 1:num_rand_samples
            sample = mvnrnd( mean, full_variance + eye(length(xstar)).*max(full_variance(:)).*0.0001, 1);
            hs = plot( xstar, sample, '-', 'Color', colorbrew(n_sample), 'Linewidth', 1); hold on;
            %
        end
    end
    
    
    h_data_adjust = plot( X, y_adjusted, 'k+', 'Linewidth', 1.5, 'Markersize', 10, 'Color', colorbrew(3)); hold on;
    
    
    %ylim( ylimits);
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
        %filename = sprintf('%s/%s-%d', introfigsdir, 'fuzzy', N );
        %myaa('publish');
        %savepng(gcf, filename);
        save2pdf([filename '.pdf'], gcf, dpi, true );
        
        % Plot separate legend
   %     if d == D
        if false
            set(h_data_orig,'Visible','off');
            set(h_data_adjust,'Visible','off');
            set(h_gp_post,'Visible','off');
            axis off
            h_legend = legend([h_data_orig, h_data_adjust, h_gp_post], ...
               'Original data', 'Adjusted', 'GP posterior', 'Location', 'Best');
            set(h_legend,'FontSize',12);
            set(h_legend,'FontName','Times');
        
            %set(gcf, 'color', 'white');
            save2pdf([filename '-legend.pdf'], gcf, dpi, true );
        end
    end
end


function col = color_spectrum(p)
    no_col = [1 1 1];
    full_col = [ 1 0 0 ];
    col = (1 - p)*no_col + p*full_col;
end

