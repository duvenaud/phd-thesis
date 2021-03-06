function plot_additive_decomp_cov( X, y, kernel_components, kernel_params, log_noise, savefigs, file_prefix, column_names )
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
dpi = 300;


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

cmax = -Inf; cmin = Inf;
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

        % Compute ranges of test sets.
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
       
        % Compute posterior covariance between these two components.
        component1_sigma_star = feval(cur_cov1{:}, cur_hyp1, x1star, X);
        component1_sigma_starstar = feval(cur_cov1{:}, cur_hyp1, x1star);
        component2_sigma_star = feval(cur_cov2{:}, cur_hyp2, x2star, X);
        component2_sigma_starstar = feval(cur_cov2{:}, cur_hyp2, x2star);
        
        covar = -component1_sigma_star / complete_sigma * component2_sigma_star';
        var1 = -component1_sigma_star / complete_sigma * component1_sigma_star';
        var2 = -component2_sigma_star / complete_sigma * component2_sigma_star';
        
        % Diagonals have an extra term.
        if i == j
            covar = covar + component1_sigma_starstar; % is equal to component2_sigma_starstar
            assert(all(component1_sigma_starstar(:) == component2_sigma_starstar(:)));
            var1 = var1 + component1_sigma_starstar;
            var2 = var2 + component2_sigma_starstar;
        end
        
        correlation = covar ./ sqrt(bsxfun(@times, diag(var1), diag(var2)'));
        
        assert(all(abs(correlation(:)) <= 1));

        % Plot posterior mean and variance.
        if savefigs
            clf; imagesc(correlation);
            caxis([-0.5 1]);
            
            % Make the plot look nice.
            set( gca, 'XTick', [] );
            set( gca, 'yTick', [] );
            set( gca, 'XTickLabel', '' );
            set( gca, 'yTickLabel', '' );
            set(gcf, 'color', 'white');
            tightfig();
            set_fig_units_cm(4,4);
            
            title1 = strrep(column_names{i}, ' ', '-');
            title2 = strrep(column_names{j}, ' ', '-');
            filename = sprintf('%s-%s', title1, title2);
            save2pdf([file_prefix, '-', filename '.pdf'], gcf, dpi, true );
            filenames{i}{j} = filename;
        else
            subplot(num_components, num_components, j + (i - 1) * num_components);
            imagesc(correlation);
            caxis([-0.5 1]);
        end
        
        % Keep track of max and min for unified colobar.
        cmax = max([cmax, max(correlation(:))]);
        cmin = min([cmin, min(correlation(:))]);
    end
end

% Print colorbar
figure;  caxis([-0.5 1]); hc = colorbar('vert'); axis off;
set(hc, 'FontName', 'Times New Roman', 'Fontsize', 12);
%tightfig();
set_fig_units_cm(10,4);
save2pdf([file_prefix, '-colorbar.pdf'], gcf, dpi, true );


if savefigs
    smalltext = '';
    if num_components > 7; smalltext = '\small'; end

    % Print a latex table.
    fprintf('\n\\renewcommand{\\tabcolsep}{1mm}');
    fprintf('\n\\def\\incpic#1{\\includegraphics[width=%1.3f\\columnwidth]{%s-#1}}', ...
            1/(num_components + 1), file_prefix);
    %fprintf('\n\\begin{tabular}{p{2mm}*{%d}{p{%1.3f\\columnwidth}}}\n & ', ...
    %        num_components, 1/(num_components + 1));
    fprintf('\n\\begin{tabular}{p{2mm}*{%d}{c}}\n & ', ...
            num_components);
    for i = 1:num_components
        % print header
        if i == 1
            for j = 1:num_components
                fprintf('%s{%s}', smalltext, column_names{j});
                if j < num_components; fprintf(' & '); else fprintf(' \\\\ \n '); end
            end
        end
        fprintf('\\rotatebox{90}{%s{%s}} & ', smalltext, column_names{i});
        for j = 1:num_components
            fprintf('\\incpic{%s}', filenames{i}{j});
            if j < num_components; fprintf(' & '); else fprintf(' \\\\ \n '); end
        end
    end
    fprintf('\\end{tabular}\n');
end

end


