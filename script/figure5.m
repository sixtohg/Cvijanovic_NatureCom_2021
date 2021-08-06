%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 5:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all;fclose all;clc
cd('/home/sixto/Documents/MLToolbox_R2013/MeteoLab/');init;
addpath('/home/sixto/Documents/Publications/Articulos/2020_Turco/lastFigure/cbrewer/cbrewer/cbrewer')
workPath='/home/sixto/Documents/Publications/Articulos/2020_Turco/';cd([workPath 'lastFigure/'])

%% log[BA] = -1.3063 + 0.1591 * Time - 0.3676 * SPI(June to October) + 0.6806 * TSMAX(April to September) + epsilon;
B1=[-1.3063 0.1591 -0.3676 0.6806];
changes_spi=-1.0:0.1:1.0;changes_tmax=-2.0:0.20:4.0;
[a1,a2]=meshgrid(changes_spi,changes_tmax);

%% CMIP5:
dir_sim_pr='/home/sixto/Documents/Publications/Articulos/2020_Turco/scripts/GB-models/CMIP5_CMIP6/CMIP5_CMIP6/pr/rcp85/';
dir_sim_tasmax='/home/sixto/Documents/Publications/Articulos/2020_Turco/data/tasmax/rcp85/';
listGCMs_SPI = dir(fullfile(dir_sim_pr, 'pr_*.mat'));
listGCMs_TX = dir(fullfile(dir_sim_tasmax, 'tasmax*_max.mat'));
indGCMs=[2 2 3 11 11 22 28 28 33 34 35 59 59 41 48 49 52 56 57 57 37];
listGCMs_SPI=listGCMs_SPI(indGCMs);
listGCMs_TX=listGCMs_TX([1 1 15 2 2 3 4 4 5 6 7 16 16 9 10 11 12 13 14 14 8]);
refTmax=repmat(NaN,20,length(listGCMs_SPI));
refSPI=repmat(NaN,20,length(listGCMs_SPI));
for ifile=1:length(listGCMs_SPI)
	load(sprintf('/home/sixto/Dropbox/ivana_ice_fires/review/data/CMIP5/spi4_3_%d_nostd_baileys.mat',ifile-1));
	load(sprintf('/home/sixto/Dropbox/ivana_ice_fires/review/data/CMIP5/tx_4_9_%d_nostd_baileys.mat',ifile-1));
    tx_4_9(find(isnan(tx_4_9.*spi4_3)))=NaN;spi4_3(find(isnan(tx_4_9.*spi4_3)))=NaN;
	refTmax(:,ifile)=tx_4_9(1:20);refSPI(:,ifile)=spi4_3(1:20);
end
refTmax=nanmean(refTmax,2);refSPI=nanmean(refSPI,2);
auxP=exp(B1(3)*refSPI+ B1(4)*refTmax);
changes_neutral=repmat(NaN,size(a1));
for i=1:size(a1,1)
	for j=1:size(a1,2)
		auxF=exp(B1(3)*(refSPI + a1(i,j))+ B1(4)*(refTmax + a2(i,j)));
		changes_neutral(i,j)=100*(nanmean(auxF)-nanmean(auxP))/nanmean(auxP);
	end
end
changes_neutral=changes_neutral';

%% CMIP5:
dir_sim_pr='/home/sixto/Documents/Publications/Articulos/2020_Turco/scripts/GB-models/CMIP5_CMIP6/CMIP5_CMIP6/pr/rcp85/';
dir_sim_tasmax='/home/sixto/Documents/Publications/Articulos/2020_Turco/data/tasmax/rcp85/';
listGCMs_SPI = dir(fullfile(dir_sim_pr, 'pr_*.mat'));
listGCMs_TX = dir(fullfile(dir_sim_tasmax, 'tasmax*_max.mat'));
indGCMs=[2 2 3 11 11 22 28 28 33 34 35 59 59 41 48 49 52 56 57 57 37];
listGCMs_SPI=listGCMs_SPI(indGCMs);
listGCMs_TX=listGCMs_TX([1 1 15 2 2 3 4 4 5 6 7 16 16 9 10 11 12 13 14 14 8]);
changes_cmip5=repmat(NaN,length(listGCMs_SPI),3);
for ifile=1:length(listGCMs_SPI)
	load(sprintf('/home/sixto/Dropbox/ivana_ice_fires/review/data/CMIP5/spi4_3_%d_nostd_baileys.mat',ifile-1));
	load(sprintf('/home/sixto/Dropbox/ivana_ice_fires/review/data/CMIP5/tx_4_9_%d_nostd_baileys.mat',ifile-1));
    tx_4_9(find(isnan(tx_4_9.*spi4_3)))=NaN;spi4_3(find(isnan(tx_4_9.*spi4_3)))=NaN;
    changes_cmip5(ifile,1)=nanmean(tx_4_9(21:40))-nanmean(tx_4_9(1:20));
    changes_cmip5(ifile,2)=nanmean(spi4_3(21:40))-nanmean(spi4_3(1:20));
	auxP=exp(B1(3)*spi4_3(1:20)+ B1(4)*tx_4_9(1:20));
	auxF=exp(B1(3)*spi4_3(21:40)+ B1(4)*tx_4_9(21:40));
	changes_cmip5(ifile,3)=100*(nanmean(auxF)-nanmean(auxP))/nanmean(auxP);
end

%% CMIP6:
dir_sim_pr='/home/sixto/Documents/Publications/Articulos/2020_Turco/scripts/GB-models/CMIP5_CMIP6/CMIP5_CMIP6/pr/ssp585/';
dir_sim_tasmax='/home/sixto/Documents/Publications/Articulos/2020_Turco/data/tasmax/ssp585/';
listGCMs_SPI = dir(fullfile(dir_sim_pr, 'pr_*.mat'));
listGCMs_TX = dir(fullfile(dir_sim_tasmax, 'tasmax*_max.mat'));
changes_cmip6=repmat(NaN,length(listGCMs_SPI),3);
for ifile=1:length(listGCMs_SPI)
	load(sprintf('/home/sixto/Dropbox/ivana_ice_fires/review/data/CMIP6/spi4_3_%d_nostd_baileys.mat',ifile-1));
	load(sprintf('/home/sixto/Dropbox/ivana_ice_fires/review/data/CMIP6/tx_4_9_%d_nostd_baileys.mat',ifile-1));
    tx_4_9(find(isnan(tx_4_9.*spi4_3)))=NaN;spi4_3(find(isnan(tx_4_9.*spi4_3)))=NaN;
    changes_cmip6(ifile,1)=nanmean(tx_4_9(21:40))-nanmean(tx_4_9(1:20));
    changes_cmip6(ifile,2)=nanmean(spi4_3(21:40))-nanmean(spi4_3(1:20));
	auxP=exp(B1(3)*spi4_3(1:20)+ B1(4)*tx_4_9(1:20));
	auxF=exp(B1(3)*spi4_3(21:40)+ B1(4)*tx_4_9(21:40));
	changes_cmip6(ifile,3)=100*(nanmean(auxF)-nanmean(auxP))/nanmean(auxP);
end

%% Low-Ice:
changes_low=repmat(NaN,12,3);
for ifile=1:12
	load(sprintf('/home/sixto/Dropbox/ivana_ice_fires/review/data/low_ice/spi4_3_%d_nostd_baileys.mat',ifile-1));
	load(sprintf('/home/sixto/Dropbox/ivana_ice_fires/review/data/low_ice/tx_4_9_%d_nostd_baileys.mat',ifile-1));
    tx_4_9(find(isnan(tx_4_9.*spi4_3)))=NaN;spi4_3(find(isnan(tx_4_9.*spi4_3)))=NaN;
    changes_low(ifile,1)=nanmean(tx_4_9(21:40))-nanmean(tx_4_9(1:20));
    changes_low(ifile,2)=nanmean(spi4_3(21:40))-nanmean(spi4_3(1:20));
	auxP=exp(B1(3)*spi4_3(1:20)+ B1(4)*tx_4_9(1:20));
	auxF=exp(B1(3)*spi4_3(21:40)+ B1(4)*tx_4_9(21:40));
	changes_low(ifile,3)=100*(nanmean(auxF)-nanmean(auxP))/nanmean(auxP);
end

%% Low-Ice EC-Earth:
changes_EC=repmat(NaN,10,3);
for ifile=1:10
	load(sprintf('/home/sixto/Dropbox/ivana_ice_fires/review/data/low_ice_ecearth/spi4_3_%d_nostd_baileys.mat',ifile-1));
	load(sprintf('/home/sixto/Dropbox/ivana_ice_fires/review/data/low_ice_ecearth/tx_4_9_%d_nostd_baileys.mat',ifile-1));
    tx_4_9(find(isnan(tx_4_9.*spi4_3)))=NaN;spi4_3(find(isnan(tx_4_9.*spi4_3)))=NaN;
    changes_EC(ifile,1)=nanmean(tx_4_9(21:40))-nanmean(tx_4_9(1:20));
    changes_EC(ifile,2)=nanmean(spi4_3(21:40))-nanmean(spi4_3(1:20));
	auxP=exp(B1(3)*spi4_3(1:20)+ B1(4)*tx_4_9(1:20));
	auxF=exp(B1(3)*spi4_3(21:40)+ B1(4)*tx_4_9(21:40));
	changes_EC(ifile,3)=100*(nanmean(auxF)-nanmean(auxP))/nanmean(auxP);
end

numcolors=11;
fall=flipud(cbrewer('div', 'PuOr', numcolors));
auxValues=linspace(-50,450,numcolors);auxAxis=(auxValues(end)-auxValues(1))/(numcolors-1);
fall(find(fall < 0))=0;fall(find(fall > 1))=1;
figure;h=surface(changes_tmax,changes_spi,changes_neutral);
axis square
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
caxis([-75 475]);colormap(fall);h=colorbar;hold on

set(h, 'YTick', [-50:50:450])
set(h,'YTickLabel',[-50:50:450])
title(h,'Changes in BA','FontWeight','bold')

set(gca,'Fontsize',12);
set(findall(gcf,'-property','FontSize'),'FontSize',12)
ylabel('Changes in SPI','Fontsize',10,'FontWeight','bold');
xlabel('Changes in TSMAX (\circC)','Fontsize',10,'FontWeight','bold');

figure,
hold on,scatter(changes_low(:,1),changes_low(:,2),100,'Marker','^','MarkerEdgeColor','black')
hold on,scatter(changes_EC(:,1),changes_EC(:,2),100,'Marker','v','MarkerEdgeColor','black')
scatter(changes_cmip5(:,1),changes_cmip5(:,2),100,'Marker','o','MarkerEdgeColor','black')
hold on,scatter(changes_cmip6(:,1),changes_cmip6(:,2),100,'Marker','s','MarkerEdgeColor','black')
caxis([-75 475]);colormap(fall);h=colorbar;hold on;%% caxis([auxValues(aux2)-10 auxValues(end)+10]);colormap(fall);h=colorbar;
legend({'Low-Ice','Low-Ice EC-Earth','CMIP5','CMIP6'},'Location','Best')

set(gca,'Fontsize',12,'xlim',[-2 4],'ylim',[-1.0 1.0]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)
ylabel('Changes in SPI','Fontsize',10,'FontWeight','bold');
xlabel('Changes in TSMAX (\circC)','Fontsize',10,'FontWeight','bold');

figure,
hold on,scatter(changes_low(:,1),changes_low(:,2),100,changes_low(:,3),'filled','Marker','^','MarkerEdgeColor','black')
hold on,scatter(changes_EC(:,1),changes_EC(:,2),100,changes_EC(:,3),'filled','Marker','^','MarkerEdgeColor','black')
scatter(changes_cmip5(:,1),changes_cmip5(:,2),100,changes_cmip5(:,3),'filled','Marker','o','MarkerEdgeColor','black')
hold on,scatter(changes_cmip6(:,1),changes_cmip6(:,2),100,changes_cmip6(:,3),'filled','Marker','s','MarkerEdgeColor','black')
caxis([-75 475]);colormap(fall);h=colorbar;hold on
legend({['Low-Ice (' num2str(median(changes_low(:,3)),'%5.1f') ')'],['Low-Ice EC-Earth (' num2str(median(changes_EC(:,3)),'%5.1f') ')'],['CMIP5 (' num2str(median(changes_cmip5(:,3)),'%5.1f') ')'],['CMIP6 (' num2str(median(changes_cmip6(:,3)),'%5.1f') ')']},'Location','Best')
%% grid on
set(h, 'YTick', [-50:50:450])
set(h,'YTickLabel',[-50:50:450])
title(h,'Changes in BA','FontWeight','bold')
set(gca,'Fontsize',12,'xlim',[-2 4],'ylim',[-1.0 1.0]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)
ylabel('Changes in SPI','Fontsize',10,'FontWeight','bold');
xlabel('Changes in TSMAX (\circC)','Fontsize',10,'FontWeight','bold');
