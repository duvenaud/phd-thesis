function p = plot_little_circles(x, y, circle, col, alpha, fast)
%AG_PLOT_LITTLE_CIRCLES Plot circular circles of relative size circle
% Returns handles to all patches plotted

if nargin < 6
    fast = false;
end

    [n, d] = size(x);
    assert( d == 1);
    assert( size(y, 2) == 1);
    
    if size(col, 2) == 1
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

    if fast
        %for k = 1:length(x)
        %    p(k) = plot(x(k), y(k), '.', 'Color', col(k,:));
        %end        
        set(0,'DefaultAxesColorOrder',jet(length(x)))
        p = plot(x, y, '.');
    else
        for k = 1:length(x)
            p(k) = patch(x(k)+mx, y(k)+my, col(k,:), 'FaceColor', col(k,:), 'FaceAlpha', alpha, 'EdgeColor', 'none');
        end
    end
end
