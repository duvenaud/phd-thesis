function p = plot_little_circles_3d(x, y, z, circle, col)
%AG_PLOT_LITTLE_CIRCLES Plot circular circles of relative size circle
% Returns handles to all patches plotted

if nargin < 6
    fast = false;
end

    [n, d] = size(x);
    assert( d == 1);
    assert( size(y, 2) == 1);
    
    if size(col, 1) == 1
        col = repmat( col, n, 1);
    end
    
    % aspect is width / height
    %fPos = get(gcf, 'Position');
    % need width, height in data values
    xl = xlim();
    yl = ylim();
    w = circle*(xl(2)-xl(1));%/fPos(3);
    h = circle*(yl(2)-yl(1));%/fPos(4);

    theta = 0:(2*pi/5):2*pi;
    mx = w*sin(theta);
    my = h*cos(theta);

    %if fast
        for k = 1:length(x)
            p(k) = plot3(x(k), y(k), z(k), '.', 'Color', col(k,:)); hold on;
        end        
    %else
    %    for k = 1:length(x)
    %        p(k) = patch(x(k)+mx, y(k)+my, col(k,:), 'FaceColor', col(k,:), 'FaceAlpha', alpha, 'EdgeColor', 'none');
    %    end
    %end
end
