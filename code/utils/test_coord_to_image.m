% test coord to color

grid_resolution = 100;
xrange = linspace( 0, 2,grid_resolution);
[x1, x2] = ndgrid( xrange, xrange );   % Make a grid
x = [ x1(:), x2(:) ];
colors = coord_to_image(x);

imagesc(reshape(colors, [grid_resolution, grid_resolution, 3]));