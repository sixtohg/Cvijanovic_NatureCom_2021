clear all
close all


%% load dati
dir_fig='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/paper/figures_temp/';
dir_out='/Users/marco/Documents/output/ivana_ice_fires/';
savename=[dir_out 'SP_changes_lowice.mat']
load(savename,'SP_changes')
savename=[dir_out 'TX_changes_lowice.mat']
load(savename,'TX_changes')
%boxplot(SP_changes)

savename=[dir_out 'SP_changes_lowice_ecearth.mat']
load(savename,'SP_changes_ecearth')
savename=[dir_out 'TX_changes_lowice_ecearth.mat']
load(savename,'TX_changes_ecearth')

savename=[dir_out 'SP_changes_rcp85.mat']
load(savename,'SP_changes_rcp85')
SP_changes_rcp85_cmip5=SP_changes_rcp85;
savename=[dir_out 'TX_changes_rcp85.mat']
load(savename,'TX_changes_rcp85')
TX_changes_rcp85_cmip5=TX_changes_rcp85;

savename=[dir_out 'SP_changes_rcp85_cmip6.mat']
load(savename,'SP_changes_rcp85')
savename=[dir_out 'TX_changes_rcp85_cmip6.mat']
load(savename,'TX_changes_rcp85')


%% SPI
dataplot=zeros(1,size(SP_changes_rcp85,2),4)*NaN;

dataplot(1,1:size(SP_changes_rcp85,2),1)=SP_changes_rcp85; 
dataplot(1,1:size(SP_changes_rcp85_cmip5,2),2)=SP_changes_rcp85_cmip5; 
dataplot(1,1:size(SP_changes_ecearth,2),3)=SP_changes_ecearth; 
dataplot(1,1:size(SP_changes,2),4)=SP_changes; 

%# centimeters units
% Y = 15 %21.0;                  %# A4 paper size
% X = 29.7;                  %# A4 paper size
% xMargin = 1;               %# left/right margins from page borders
% yMargin = 1;               %# bottom/top margins from page borders
% xSize = X - 2*xMargin;     %# figure size on paper (widht & hieght)
% ySize = Y - 2*yMargin;     %# figure size on paper (widht & hieght)

fntsz=12;

c=gray(3);
c=c(2,:);

disp([prctile(dataplot(1,:,1),[2.5 5 25 50 75 95 97.5])])
disp([prctile(dataplot(1,:,2),[2.5 5 25 50 75 95 97.5])])
disp([prctile(dataplot(1,:,3),[2.5 5 25 50 75 95 97.5])])
disp([prctile(dataplot(1,:,4),[2.5 5 25 50 75 95 97.5])])



hFig=figure;aboxplot3(dataplot,'colormap',c,'labels',{'CMIP6','CMIP5','low-ice EC-EARTH','low-ice PSIPPS'})
xticklabel_rotate([],45,[]);
%ymax=1.1;
%ylim([-0.65 ymax])

%RCMs+MLR uncertainties RCMs uncertainties
ylabel('Ensemble mean differences','FontSize',fntsz);
xlabel('Scenario','FontSize',fntsz);
legend({'SPI'},'Location','northoutside','Orientation','horizontal','FontSize',fntsz)

gridxy([],[-1:0.5:1],'Color','k','Linestyle',':');
gridxy([],[0],'Color','k','Linestyle','-');
set(gca,'FontSize',fntsz)

%# figure size displayed on screen (50% scaled, but same aspect ratio)
%set(hFig, 'Units','centimeters', 'Position',[0 0 xSize ySize]/2)
%movegui(hFig, 'center')

%# figure size printed on paper
% set(hFig, 'PaperUnits','centimeters')
% set(hFig, 'PaperSize',[X Y])
% set(hFig, 'PaperPosition',[xMargin yMargin xSize ySize])
% set(hFig, 'PaperOrientation','portrait')

nomeout=[dir_fig,'figure2b_temp.eps'];
print( gcf, '-depsc2', nomeout ,'-painters')


%% TX

dataplot=zeros(1,size(TX_changes_rcp85,2),4)*NaN;

dataplot(1,1:size(TX_changes_rcp85,2),1)=TX_changes_rcp85; 
dataplot(1,1:size(TX_changes_rcp85_cmip5,2),2)=TX_changes_rcp85_cmip5; 
dataplot(1,1:size(TX_changes_ecearth,2),3)=TX_changes_ecearth; 
dataplot(1,1:size(TX_changes,2),4)=TX_changes; 


%# centimeters units
% Y = 15 %21.0;                  %# A4 paper size
% X = 29.7;                  %# A4 paper size
% xMargin = 1;               %# left/right margins from page borders
% yMargin = 1;               %# bottom/top margins from page borders
% xSize = X - 2*xMargin;     %# figure size on paper (widht & hieght)
% ySize = Y - 2*yMargin;     %# figure size on paper (widht & hieght)

fntsz=12;

c=gray(3);
c=c(2,:);

disp([prctile(dataplot(1,:,1),[2.5 5 25 50 75 95 97.5])])
disp([prctile(dataplot(1,:,2),[2.5 5 25 50 75 95 97.5])])
disp([prctile(dataplot(1,:,3),[2.5 5 25 50 75 95 97.5])])
disp([prctile(dataplot(1,:,4),[2.5 5 25 50 75 95 97.5])])

hFig=figure;aboxplot3(dataplot,'colormap',c,'labels',{'CMIP6','CMIP5','low-ice EC-EARTH','low-ice PSIPPS'})
xticklabel_rotate([],45,[]);
set(gca,'FontSize',fntsz)

% ymax=3;
% ylim([-1 ymax])




%RCMs+MLR uncertainties RCMs uncertainties
ylabel('Ensemble mean differences (\circC)','FontSize',fntsz);
xlabel('Scenario','FontSize',fntsz);
legend({'TMAX'},'Location','northoutside','Orientation','horizontal','FontSize',fntsz)

%gridxy([],[-1:0.5:ymax],'Color','k','Linestyle',':');
gridxy([],[0],'Color','k','Linestyle','-');
set(gca,'FontSize',fntsz)

%# figure size displayed on screen (50% scaled, but same aspect ratio)
%set(hFig, 'Units','centimeters', 'Position',[0 0 xSize ySize]/2)
%movegui(hFig, 'center')

%# figure size printed on paper
% set(hFig, 'PaperUnits','centimeters')
% set(hFig, 'PaperSize',[X Y])
% set(hFig, 'PaperPosition',[xMargin yMargin xSize ySize])
% set(hFig, 'PaperOrientation','portrait')

nomeout=[dir_fig,'figure2a_temp.eps'];
print( gcf, '-depsc2', nomeout ,'-painters')
