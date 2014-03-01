function set_fig_units_cm( width, height )

set(gcf,'units','centimeters');
pos = get(gcf,'position');
set(gcf,'position',[pos(1:2),width .* 1.4 ,height .* 1.4]);
