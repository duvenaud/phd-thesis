% Plot the density in the observed space.
%
% David and Tomo

function draw_warped_density( X, Y, log_hypers, ...
                              circle_size, circle_alpha,N_points)

figure(5432100); clf;

[N,Q] = size(X);
[N,D] = size(Y);

if nargin < 4
    circle_size = 0.05;
    circle_alpha = 0.2;
    N_points = 1000;
end

latent_draws = mvnrnd( zeros(1, Q), eye(Q), N_points);

% Compute Gram gram matrix.
hyp(1) = -log_hypers.gamma/2;
hyp(2) = log_hypers.alpha/2;
K = covSEiso(hyp, X) + eye(N)*max(exp(log_hypers.betainv), 1e-3);  % Add a little noise.

% Compute conditional posterior.
crosscov = covSEiso(hyp, latent_draws, X);
post_mean = crosscov*(K\Y);
prior_var = covSEiso(hyp, latent_draws, 'diag');
post_var = prior_var - sum(bsxfun(@times, crosscov/K, crosscov), 2);

samples = post_mean + randn(N_points, D) .* repmat(post_var, 1, D);


plot_little_circles(samples(:,1), samples(:,2), ...
                    circle_size, colorbrew(3), circle_alpha );
hold on;

% Draw the original data and current assigments
% markers = {'x', 'o', 's', 'd', '.', '>', '<', '^'};
plot(Y(:, 1), Y(:, 2), 'o', 'MarkerEdgeColor', 'k');
hold on;

title('Observed space');
end


