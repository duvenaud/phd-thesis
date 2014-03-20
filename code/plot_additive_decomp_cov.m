function plot_additive_decomp_cov( X, y, kernel_components, kernel_params, log_noise, show_samples, savefigs, fileprefix )
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



figure; clf;
%figure(1); clf;
%filename = sptrintf( '%s-complete', fileprefix );
%nice_oned_plot( X, y, xstar, complete_mean, complete_var, true, filename, savefigs)

left_extend = 0.05;
right_extend = 0.05;

% Next, show the posterior of each component, one at a time.
for i = 1:num_components
    for j = 1:num_components

        cur_cov1 = kernel_components{i};
        cur_hyp1 = kernel_params{i};    
        cur_cov2 = kernel_components{j};
        cur_hyp2 = kernel_params{j};    

        % Figure out which dimension this kernel applies to.
        assert( strcmp( cur_cov1{1}, 'covMask' ) );
        d1 = cur_cov1{2}{1};
        assert( strcmp( cur_cov2{1}, 'covMask' ) );
        d2 = cur_cov2{2}{1};


        x1_left = min(X(:,d1)) - (max(X(:,d1)) - min(X(:,d1)))*left_extend;
        x1_right = max(X(:,d1)) + (max(X(:,d1)) - min(X(:,d1)))*right_extend;

        x2_left = min(X(:,d2)) - (max(X(:,d2)) - min(X(:,d2)))*left_extend;
        x2_right = max(X(:,d2)) + (max(X(:,d2)) - min(X(:,d2)))*right_extend;

        x1range = linspace(x1_left, x1_right, nstar)';
        x1star = NaN(nstar,D);
        x1star(:, d1) = x1range;

        x2range = linspace(x2_left, x2_right, nstar)';
        x2star = NaN(nstar,D);
        x2star(:, d2) = x2range;    




        % Compute Gram matrices of just this component.
        component1_sigma = feval(cur_cov1{:}, cur_hyp1, X);
        component1_sigma_star = feval(cur_cov1{:}, cur_hyp1, x1star, X);
        component1_sigma_starstar = feval(cur_cov1{:}, cur_hyp1, x1star);

        component2_sigma = feval(cur_cov2{:}, cur_hyp2, X);
        component2_sigma_star = feval(cur_cov2{:}, cur_hyp2, x2star, X);
        component2_sigma_starstar = feval(cur_cov2{:}, cur_hyp2, x2star);    

        covar = component1_sigma_star / complete_sigma * component2_sigma_star';
        
        % diagonals have an extra term
        if i == j
            covar = covar + component1_sigma_starstar;
        end


        % Plot posterior mean and variance.
        subplot(num_components, num_components, j + (i - 1) * num_components);
        imagesc(covar);
    end
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

