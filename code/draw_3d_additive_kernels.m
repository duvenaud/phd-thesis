function draw_3d_additive_kernels
% Make figures of additive kernels in 3 dimensions.
%
% David Duvenaud
% March 2014
% =================

addpath(genpath([pwd '/gpml/']))
addpath('utils/');

clear all
close all

dpi = 200;

% X,Y,Z iz the meshgrid and V is my function evaluated at each meshpoint
range = -17:.3:17;
[X,Y,Z] = meshgrid(range, range, range );
xstar = [X(:) Y(:) Z(:)];

a = .9;

first_order_variance = 1.07;
second_order_variance = 25;
third_order_variance = 1000009;

set_fig_units_cm(5,5);
covfunc = {'covADD',{[1,2,3],'covSEiso'}};
hyp.cov = log([1,1,1,1,1,1,first_order_variance, second_order_variance, third_order_variance]);
V = feval(covfunc{:}, hyp.cov, xstar, [0,0,0]);
V = reshape(V, length(range),length(range),length(range));
figure(1);  draw_isosurface( X, Y, Z, V, a, range);
nice_figure_save('3d_add_kernel_321');

covfunc = {'covADD',{[2,3],'covSEiso'}};
hyp.cov = log([1,1,1,1,1,1, second_order_variance, third_order_variance]);
V = feval(covfunc{:}, hyp.cov, xstar, [0,0,0]);
V = reshape(V, length(range),length(range),length(range));
figure(2); draw_isosurface( X, Y, Z, V, a, range);
nice_figure_save('3d_add_kernel_32');

covfunc = {'covADD',{[3],'covSEiso'}};
hyp.cov = log([1,1,1,1,1,1, third_order_variance]);
V = feval(covfunc{:}, hyp.cov, xstar, [0,0,0]);
V = reshape(V, length(range),length(range),length(range));
figure(3); draw_isosurface( X, Y, Z, V, a, range);
nice_figure_save('3d_add_kernel_3');

covfunc = {'covADD',{[2],'covSEiso'}};
hyp.cov = log([1,1,1,1,1,1, 25]);
V = feval(covfunc{:}, hyp.cov, xstar, [0,0,0]);
V = reshape(V, length(range),length(range),length(range));
figure(4); draw_isosurface( X, Y, Z, V, a, range);
nice_figure_save('3d_add_kernel_2');

covfunc = {'covADD',{[1],'covSEiso'}};
hyp.cov = log([1,1,1,1,1,1,first_order_variance]);
V = feval(covfunc{:}, hyp.cov, xstar, [0,0,0]);
V = reshape(V, length(range),length(range),length(range));
figure(5); draw_isosurface( X, Y, Z, V, a, range);
nice_figure_save('3d_add_kernel_1');


end

function draw_isosurface( X, Y, Z, V, a, range)
    p = patch(isosurface(X,Y,Z,V,a)); % isosurfaces at max(V)/a
    isonormals(X,Y,Z,V,p); % plot the surfaces
    set(p,'FaceColor','red','EdgeColor','none'); % set colors
    camlight; camlight(-80,-10); lighting gouraud; 
    %alpha(.1); % set the transparency for the isosurfaces
    view(3); daspect([1 1 1]); %axis tight;
    axis( [ min(range), max(range),min(range), max(range),min(range), max(range)]);
    %axis off;
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
    savepng(gcf, ['../figures/additive/3d-kernel/', filename]);    
    %filename_eps = ['../figures/additive/3d-kernel/', filename, '.eps']
    %filename_pdf = ['../figures/additive/3d-kernel/', filename, '.pdf']
    %print -depsc2 filename_eps
    %eps2pdf( filename_eps, filename_pdf, true);
end