function plot_satnumlog(satnumlog)

figure(30)
clf;
scatter(satnumlog(:,4), satnumlog(:,3), 50,satnumlog(:,6),'filled');


figure(31)
clf;
plot(satnumlog(:,1), satnumlog(:,6), '-', 'linewidth',2);