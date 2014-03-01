% test coord to color

n = 1000;
x = randn( n, 2);
colors = coord_to_color2(x);

%close all;
figure;
for i = 1:n
    plot(x(i,1), x(i,2), '.', 'color', colors(i,:)); hold on;
end