function plot_additive_decomp( X, y, kernel_components, kernel_params, log_noise, show_samples, savefigs, fileprefix )
%
% Decomposes an additive GP model into component parts
%
% X,y: The original data
% kernel_components: a cell array of covariance functions, to be added
%                    together.  This code can figure out which dimension
%                    each kernel applies to.
% kernel_params: a cell array of kernel hyperparameters.
%
% David Duvenaud
% March 2014

addpath(genpath( 'utils' ));
addpath(genpath( 'gpml' ));

if nargin < 7; savefigs = false; end
if nargin < 6; show_samples = true; end

nstar = 400;




[N,D] = size(X);
num_components = length(kernel_components);
assert(num_components == length(kernel_params));

noise_var = exp(2*log_noise);

% First, compute and plot the entire model.
% concatenate all additive components
complete_cov = { 'covSum', kernel_components };
complete_hypers = unwrap( kernel_params );

noise_cov = eye(length(y)).*noise_var;
complete_sigma = feval(complete_cov{:}, complete_hypers, X) + noise_cov;

%complete_sigmastar = feval(complete_cov{:}, complete_hypers, X, xrange);
%complete_sigmastarstart = feval(complete_cov{:}, complete_hypers, xrange);

% First, plot the data
%complete_mean = complete_sigmastar' / complete_sigma * y;
%complete_var = complete_sigmastarstart - complete_sigmastar' / complete_sigma * complete_sigmastar;

% Detect categorical dimensions.
categorical = zeros(1,D);
for d = 1:D
    num_unique = length(unique(X(:,d)));
    if num_unique < 10
        categorical(d) = true;
        fprintf('Dimension %d is categorical with %d different values\n', ...
                d, num_unique);
    end
end




%figure(1); clf;
%filename = sptrintf( '%s-complete', fileprefix );
%nice_oned_plot( X, y, xstar, complete_mean, complete_var, true, filename, savefigs)



% Next, show the posterior of each component, one at a time.
for i = 1:num_components
    
    cur_cov = kernel_components{i};
    cur_hyp = kernel_params{i};    
    
    % Figure out which dimension this kernel applies to.
    assert( strcmp( cur_cov{1}, 'covMask' ) );
    d = cur_cov{2}{1};
    
    for zoom = [true, false]
        if zoom
            left_extend = -0.5;
            right_extend = -0.49;
        else
            left_extend = 0.05;
            right_extend = 0.2;
        end
        x_left = min(X(:,d)) - (max(X(:,d)) - min(X(:,d)))*left_extend;
        x_right = max(X(:,d)) + (max(X(:,d)) - min(X(:,d)))*right_extend;
        xrange = linspace(x_left, x_right, nstar)';
        xstar = NaN(nstar,D);
        xstar(:, d) = xrange;



        % Compute Gram matrices of just this component.
        component_sigma = feval(cur_cov{:}, cur_hyp, X);
        component_sigma_star = feval(cur_cov{:}, cur_hyp, xstar, X);
        component_sigma_starstar = feval(cur_cov{:}, cur_hyp, xstar);

        % Compute posterior of just this component.
        component_mean = component_sigma_star / complete_sigma * y;
        component_var = component_sigma_starstar - component_sigma_star / complete_sigma * component_sigma_star';

        data_mean = component_sigma / complete_sigma * y;
        %data_mean = y - (complete_sigma - component_sigma) / complete_sigma * y;


        % Plot posterior mean and variance.
        figure(i); clf;
        filename = sprintf( '%s-component-%d', fileprefix, i );
        if zoom
            filename = [filename, '-zoom'];
        end
        nice_oned_plot( X(:,d), y, data_mean, xrange, component_mean, component_var, show_samples, savefigs, filename)
    end
end

complete_mean = (complete_sigma - noise_cov) / complete_sigma * y;
rs = 1 - var(complete_mean - y)/ var(y)



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
    
    h_data_orig = plot( X, y, 'k.', 'Linewidth', 1.5, 'Markersize', 10, 'Color', colorbrew(2)); hold on;

    if show_samples
        % Plot posterior samples
        
        
        seed=0;   % fixing the seed of the random generators
        randn('state',seed);
        rand('state',seed);

        for n_sample = 1:num_rand_samples
            L = chol(full_variance + eye(length(xstar)).*max(full_variance(:)).*0.0001);
            sample = mean + L'*randn(length(xstar),1);
            %sample = mvnrnd( mean, full_variance + eye(length(xstar)).*max(full_variance(:)).*0.0001, 1);
            hs = plot( xstar, sample, '-', 'Color', colorbrew(n_sample), 'Linewidth', 1); hold on;
            %
        end
    end
    
    
    h_data_adjust = plot( X, y_adjusted, 'k+', 'Linewidth', 1.5, 'Markersize', 5, 'Color', colorbrew(3)); hold on;
    
    
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
    set_fig_units_cm(11,8);

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

