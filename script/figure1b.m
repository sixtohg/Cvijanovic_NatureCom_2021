
%% frap
clear all
close all

%%

years=1951:2020;
dir_fire='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/';
dir_out='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/paper/figures_temp/';
domain='baileys';

%% frap 
savename = [dir_fire,'fire_1950_2020_',domain,'_mjjas_gt300.mat'];
%savename = [dir_fire,'fire_1950_2020_',domain,'_mjjaso.mat'];
load(savename,'FIRE')
anni_fire=1950:2020;
[C,IA,IB] = intersect(years,anni_fire);
FIRE=FIRE(IB);

%yor=detrend(log(FIRE));




%% califonia nclimgrid
filename = [dir_fire,'spi4_3_nclimgrid.mat'];
load(filename);
filename = [dir_fire,'tx_4_9_nclimgrid.mat'];
load(filename);

spi4_3=scale_1981_2000(spi4_3,years);
tx_4_9=scale_1981_2000(tx_4_9,years);

%% scatter plot 1


BAS_plot=(FIRE); %*1000;

%%

figure1 = figure('PaperSize',[20.98 29.68]);
% Create axes
axes1 = axes('Parent',figure1);
hold('all');



scatter1 =scatter(spi4_3,tx_4_9,BAS_plot,years,'MarkerFaceColor','flat','MarkerEdgeColor',[0 0 0])
h=colorbar
cmap = (autumn(length(years)));
colormap(cmap)
title(h,'Year','FontSize',18)
xlabel('SPI','FontSize',12)
ylabel('TSMAX','FontSize',12)
scatter1.MarkerFaceAlpha = 0.8;
xlim([-3 3])
ylim([-3 3.9])


hold on

gridxy(0,'color','k')
line([-3,3],[0,0],'color','k')
%line([-3 0],[3 0],'color','k')


%ylim([7.5 14.5])
xlabel('SPI','FontSize',12)
ylabel('TSMAX','FontSize',12)
set(gca,'FontSize',12)
hold on

scatter(1,2.6,5000,'filled','black')
text(0.7,3.6,[num2str(5000) ' km^2'])
scatter(1.9,2.6,500,'filled','black')
text(1.6,3.6,[num2str(500) ' km^2'])
scatter(2.7,2.6,100,'filled','black')
text(2.4,3.6,[num2str(100) ' km^2'])

file=[dir_out,'figure1b_temp.eps']
print( gcf, '-depsc2', file ,'-painters')
set(gcf,'PaperType','A4')

print( gcf, '-dpdf', [file(1:end-4),'.pdf'] ,'-painters')