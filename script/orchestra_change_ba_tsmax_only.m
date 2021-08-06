clear all,clc, close all
%% cd('/home/sixto/Documents/MLToolbox_R2013/MeteoLab/'),init
workPath=[pwd '/../'];cd([workPath 'script/'])
dir_out='/Users/marco/Documents/output/ivana_ice_fires/';
domain='baileys';
dir_fire=[workPath 'data/'];
dir_fig=[workPath 'paper/figures_temp/']; 

version_name='tasmax';


%%

years=(1951:2020)';
years_sim=[1981:2000, 2031:2050]';

alpha = 0.05;
NB=1000;

%%
savename = [dir_fire,'fire_1950_2020_',domain,'_mjjas_gt300.mat'];
%savename = [dir_fire,'fire_1950_2020_',domain,'_mjjaso.mat'];
load(savename,'FIRE')
anni_fire=1950:2020;
[C,IA,IB] = intersect(years,anni_fire);
FIRE=FIRE(IB)/1000;


yor=(log(FIRE))';
figure; plot(yor)


%% califonia nclimgrid
filename = [dir_fire,'spi4_3_nclimgrid.mat'];
load(filename);
filename = [dir_fire,'tx_4_9_nclimgrid.mat'];
load(filename);


SP=scale_1981_2000(spi4_3,years);
TX=scale_1981_2000(tx_4_9,years);


[a,b]=corrcoef(SP,TX)

%% regression climate

T=years';

Xorig1 = [ones(size(SP,1),1) SP TX years (years.*years) ]; %osservati

[B1,BINT,R,RINT,STATS1] = regress((yor),Xorig1,alpha);
STATS1


Y1 = (Xorig1*B1);
[a,b]=corrcoef((yor),Y1)



figure;
plot(years,(yor),'-ko','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
hold on
plot(years,(Y1),'k:o','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',6)
legend('Observed','Modelled','Location','Best')
xlabel('Years') 

bootdat = [(yor),Xorig1];
bootb1 = bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat);
bootb1 = sort(bootb1);
bootCI1 = [bootb1(25,:); bootb1(975,:)]  %Report 5% confidence intervals
B1

%% CMIP6

dir_sim_pr=[workPath 'data/CMIP5_CMIP6/pr/ssp585/'];
dir_sim_tasmax=[workPath 'data/CMIP5_CMIP6/tasmax/ssp585/'];


listGCMs_SPI = dir(fullfile(dir_sim_pr, 'pr_*_spi4.mat'));
listGCMs_TX = dir(fullfile(dir_sim_tasmax, 'tasmax*.mat'));

k=0;
for ifile=1:length(listGCMs_SPI)
    
      load([dir_sim_pr listGCMs_SPI(ifile).name])
      load([dir_sim_tasmax listGCMs_TX(ifile).name])

    pred = B1(1)  + B1(3) * scale_1981_2000(tx_4_9,years_sim) + B1(4) * mean(years) + B1(5) * mean(years.*years);
    ctrl_ens=nanmean(exp(pred(1:20)));
    ice_ens=nanmean(exp(pred(21:40)));
    changes(ifile)=100*(ice_ens-ctrl_ens)/ctrl_ens;
    
    for ib=1:1000
        k=k+1;
        pred = bootb1(ib,1)  + bootb1(ib,3) * scale_1981_2000(tx_4_9,years_sim) + bootb1(ib,4) * mean(years') + bootb1(ib,5) * mean(years.*years)';
        ctrl_ens=nanmean(exp(pred(1:20)));
        ice_ens=nanmean(exp(pred(21:40)));
        changes_boot(k)=100*(ice_ens-ctrl_ens)/ctrl_ens;
    end
    
    
end

%  save
savename=[dir_out 'change_BA_rcp85_cmip6_' version_name '.mat'];
save(savename,'changes')
change_ba4=changes; clear changes
savename=[dir_out 'change_BA_boot_rcp85_cmip6_' version_name '.mat']
save(savename,'changes_boot')
change_ba4_boot=changes_boot; clear changes_boot

%% CMIP5

dir_sim_pr=[workPath 'data/CMIP5_CMIP6/pr/rcp85/'];
dir_sim_tasmax=[workPath 'data/CMIP5_CMIP6/tasmax/rcp85/'];
listGCMs={'ACCESS1-3','ACCESS1-3','BCC-CSM1-1','CNRM-CM5','CNRM-CM5','CanESM2','EC-EARTH','EC-EARTH','FGOALS-g2','GFDL-ESM2G', ...
    'GFDL-ESM2M','inmcm4','inmcm4','IPSL-CM5A-LR','MIROC-ESM','MIROC5','MPI-ESM-LR','MRI-CGCM3','NorESM1-M','NorESM1-M','HadGEM2-ES'};


listGCMs_SPI = dir(fullfile(dir_sim_pr, 'pr_*_spi4.mat'));
listGCMs_TX = dir(fullfile(dir_sim_tasmax, 'tasmax*.mat'));

k=0;
for ifile=1:length(listGCMs)
    
      load([dir_sim_pr,'spi4_3_',num2str(ifile-1),'_nostd_',domain,'.mat'])
      load([dir_sim_tasmax,'tx_4_9_',num2str(ifile-1),'_nostd_',domain,'.mat'])
      

   pred = B1(1)  + B1(3) * scale_1981_2000(tx_4_9,years_sim) + B1(4) * mean(years) + B1(5) * mean(years.*years);
    
    ctrl_ens=nanmean(exp(pred(1:20)));
    ice_ens=nanmean(exp(pred(21:40)));
    
    changes(ifile)=100*(ice_ens-ctrl_ens)/ctrl_ens;
    
    
    for ib=1:1000
        k=k+1;
        pred = bootb1(ib,1) + bootb1(ib,3) * scale_1981_2000(tx_4_9,years_sim) + bootb1(ib,4) * mean(years') + bootb1(ib,5) * mean(years.*years)';
        ctrl_ens=nanmean(exp(pred(1:20)));
        ice_ens=nanmean(exp(pred(21:40)));
        changes_boot(k)=100*(ice_ens-ctrl_ens)/ctrl_ens;
    end
    
    
end


%  save
savename=[dir_out 'change_BA_rcp85_' version_name '.mat']
save(savename,'changes')
change_ba3=changes; clear changes
savename=[dir_out 'change_BA_boot_rcp85_' version_name '.mat']
save(savename,'changes_boot')
change_ba3_boot=changes_boot; clear changes_boot

%% low-ice

dir_sim=[workPath 'data/low_ice/'];


k=0;
for ifile=1:12
    
      load([dir_sim 'spi4_3_' num2str(ifile-1) '_nostd_' char(domain) '.mat'])
      
      load([dir_sim  'tx_4_9_' num2str(ifile-1) '_nostd_' char(domain) '.mat'])

      
    pred = B1(1)  + B1(3) * scale_1981_2000(tx_4_9,years_sim) + B1(4) * mean(years) + B1(5) * mean(years.*years);
    
    ctrl_ens=nanmean(exp(pred(1:20)));
    ice_ens=nanmean(exp(pred(21:40)));
    
    changes(ifile)=100*(ice_ens-ctrl_ens)/ctrl_ens;
    
    
    for ib=1:1000
        k=k+1;
        pred = bootb1(ib,1)  + bootb1(ib,3) * scale_1981_2000(tx_4_9,years_sim) + bootb1(ib,4) * mean(years') + bootb1(ib,5) * mean(years.*years)';
        ctrl_ens=nanmean(exp(pred(1:20)));
        ice_ens=nanmean(exp(pred(21:40)));
        changes_boot(k)=100*(ice_ens-ctrl_ens)/ctrl_ens;
    end
    
    
end


%% low-ice +2
k=0;
for ifile=1:12
    
  
    load([dir_sim 'spi4_3_' num2str(ifile-1) '_nostd_' char(domain) '.mat'])
      
      load([dir_sim  'tx_4_9_' num2str(ifile-1) '_nostd_' char(domain) '.mat'])
      tx_4_9(21:end)=tx_4_9(21:end)+2;
      
     pred = B1(1) + B1(3) * scale_1981_2000(tx_4_9,years_sim) + B1(4) * mean(years) + B1(5) * mean(years.*years);
    
    ctrl_ens=nanmean(exp(pred(1:20)));
    ice_ens=nanmean(exp(pred(21:40)));
    
    changes_2(ifile)=100*(ice_ens-ctrl_ens)/ctrl_ens;
    
    
    for ib=1:1000
        k=k+1;
        pred = bootb1(ib,1)  + bootb1(ib,3) * scale_1981_2000(tx_4_9,years_sim) + bootb1(ib,4) * mean(years') + bootb1(ib,5) * mean(years.*years)';
        ctrl_ens=nanmean(exp(pred(1:20)));
        ice_ens=nanmean(exp(pred(21:40)));
        changes_2_boot(k)=100*(ice_ens-ctrl_ens)/ctrl_ens;
    end
    
    
end

%  save
savename=[dir_out 'change_BA_lowice_' version_name '.mat']
save(savename,'changes')
savename=[dir_out 'change_BA_boot_lowice_' version_name '.mat']
save(savename,'changes_boot')

savename=[dir_out 'change_2_BA_lowice_' version_name '.mat']
save(savename,'changes_2')
savename=[dir_out 'change_2_BA_boot_lowice_' version_name '.mat']
save(savename,'changes_2_boot')


%% low-ice ec-earth

dir_sim=[workPath 'data/low_ice_ecearth/'];

%% low-ice
SP=zeros(40,12)*NaN;
TX=zeros(40,12)*NaN;
k=0;
for ifile=1:10
    
      load([dir_sim 'spi4_3_' num2str(ifile-1) '_nostd_' char(domain) '.mat'])
      load([dir_sim  'tx_4_9_' num2str(ifile-1) '_nostd_' char(domain) '.mat'])
      
        pred = B1(1) + B1(3) * scale_1981_2000(tx_4_9,years_sim) + B1(4) * mean(years) + B1(5) * mean(years.*years);
    
    ctrl_ens=nanmean(exp(pred(1:20)));
    ice_ens=nanmean(exp(pred(21:40)));
    
    changes(ifile)=100*(ice_ens-ctrl_ens)/ctrl_ens;
    
    
    for ib=1:1000
        k=k+1;
        pred = bootb1(ib,1)  + bootb1(ib,3) * scale_1981_2000(tx_4_9,years_sim) + bootb1(ib,4) * mean(years') + bootb1(ib,5) * mean(years.*years)';
        ctrl_ens=nanmean(exp(pred(1:20)));
        ice_ens=nanmean(exp(pred(21:40)));
        changes_boot(k)=100*(ice_ens-ctrl_ens)/ctrl_ens;
    end
    
    
end


%% low-ice +2
k=0;
for ifile=1:10
    
  
    load([dir_sim 'spi4_3_' num2str(ifile-1) '_nostd_' char(domain) '.mat'])
      
      load([dir_sim  'tx_4_9_' num2str(ifile-1) '_nostd_' char(domain) '.mat'])
      tx_4_9(21:end)=tx_4_9(21:end)+2;
      

    pred = B1(1) + B1(3) * scale_1981_2000(tx_4_9,years_sim) + B1(4) * mean(years) + B1(5) * mean(years.*years);
    
    ctrl_ens=nanmean(exp(pred(1:20)));
    ice_ens=nanmean(exp(pred(21:40)));
    
    changes_2(ifile)=100*(ice_ens-ctrl_ens)/ctrl_ens;
    
    
    for ib=1:1000
        k=k+1;
        %pred = bootb1(ib,1) + bootb1(ib,2) * spi4_3 + bootb1(ib,3) * tx_4_9;
        pred = bootb1(ib,1)  + bootb1(ib,3) * scale_1981_2000(tx_4_9,years_sim) + bootb1(ib,4) * mean(years') + bootb1(ib,5) * mean(years.*years)';
        ctrl_ens=nanmean(exp(pred(1:20)));
        ice_ens=nanmean(exp(pred(21:40)));
        changes_2_boot(k)=100*(ice_ens-ctrl_ens)/ctrl_ens;
    end
    
    
end

%  save
savename=[dir_out 'change_BA_lowice_ecearth_' version_name '.mat']
save(savename,'changes')
savename=[dir_out 'change_BA_boot_lowice_ecearth_' version_name '.mat']
save(savename,'changes_boot')

savename=[dir_out 'change_2_BA_lowice_ecearth_' version_name '.mat']
save(savename,'changes_2')
savename=[dir_out 'change_2_BA_boot_lowice_ecearth_' version_name '.mat']
save(savename,'changes_2_boot')

%% plot


dir_fig='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/paper/figures_temp/';
%% load dati

dir_out='/Users/marco/Documents/output/ivana_ice_fires/';
savename=[dir_out 'change_BA_lowice_' version_name '.mat']
load(savename,'changes')
change_ba1=changes; clear changes
savename=[dir_out 'change_BA_boot_lowice_' version_name '.mat']
load(savename,'changes_boot')
change_ba1_boot=changes_boot; clear changes_boot

savename=[dir_out 'change_2_BA_lowice_' version_name '.mat']
load(savename,'changes_2')
change_ba2=changes_2; clear changes_2
savename=[dir_out 'change_2_BA_boot_lowice_' version_name '.mat']
load(savename,'changes_2_boot')
change_ba2_boot=changes_2_boot; clear changes_2_boot

savename=[dir_out 'change_BA_lowice_ecearth_' version_name '.mat']
load(savename,'changes')
change_ba3=changes; clear changes
savename=[dir_out 'change_BA_boot_lowice_ecearth_' version_name '.mat']
load(savename,'changes_boot')
change_ba3_boot=changes_boot; clear changes_boot

savename=[dir_out 'change_2_BA_lowice_ecearth_' version_name '.mat']
load(savename,'changes_2')
change_ba4=changes_2; clear changes_2
savename=[dir_out 'change_2_BA_boot_lowice_ecearth_' version_name '.mat']
load(savename,'changes_2_boot')
change_ba4_boot=changes_2_boot; clear changes_2_boot

savename=[dir_out 'change_BA_rcp85_' version_name '.mat']
load(savename,'changes')
change_ba5=changes; clear changes
savename=[dir_out 'change_BA_boot_rcp85_' version_name '.mat']
load(savename,'changes_boot')
change_ba5_boot=changes_boot; clear changes_boot

savename=[dir_out 'change_BA_rcp85_cmip6_' version_name '.mat']
load(savename,'changes')
change_ba6=changes; clear changes
savename=[dir_out 'change_BA_boot_rcp85_cmip6_' version_name '.mat']
load(savename,'changes_boot')
change_ba6_boot=changes_boot; clear changes_boot

%%
dataplot=zeros(3,size(change_ba6_boot,2),6)*NaN;

dataplot(1,1:size(change_ba6,2),1)=change_ba6; 
dataplot(2,1:size(change_ba6_boot,2),1)=change_ba6_boot; 
    
dataplot(1,1:size(change_ba5,2),2)=change_ba5; 
dataplot(2,1:size(change_ba5_boot,2),2)=change_ba5_boot; 

dataplot(1,1:size(change_ba4,2),3)=change_ba4; 
dataplot(2,1:size(change_ba4_boot,2),3)=change_ba4_boot; 

dataplot(1,1:size(change_ba2,2),4)=change_ba2; 
dataplot(2,1:size(change_ba2_boot,2),4)=change_ba2_boot; 


dataplot(1,1:size(change_ba3,2),5)=change_ba3; 
dataplot(2,1:size(change_ba3_boot,2),5)=change_ba3_boot; 

    
dataplot(1,1:size(change_ba1,2),6)=change_ba1; 
dataplot(2,1:size(change_ba1_boot,2),6)=change_ba1_boot; 


disp([prctile(dataplot(1,:,1),[2.5 5 25 50 75 95 97.5])])
disp([prctile(dataplot(2,:,1),[2.5 5 25 50 75 95 97.5])])

disp([prctile(dataplot(1,:,2),[2.5 5 25 50 75 95 97.5])])
disp([prctile(dataplot(2,:,2),[2.5 5 25 50 75 95 97.5])])

disp([prctile(dataplot(1,:,3),[2.5 5 25 50 75 95 97.5])])
disp([prctile(dataplot(2,:,3),[2.5 5 25 50 75 95 97.5])])

disp([prctile(dataplot(1,:,4),[2.5 5 25 50 75 95 97.5])])
disp([prctile(dataplot(2,:,4),[2.5 5 25 50 75 95 97.5])])

disp([prctile(dataplot(1,:,5),[2.5 5 25 50 75 95 97.5])])
disp([prctile(dataplot(2,:,5),[2.5 5 25 50 75 95 97.5])])

disp([prctile(dataplot(1,:,6),[2.5 5 25 50 75 95 97.5])])
disp([prctile(dataplot(2,:,6),[2.5 5 25 50 75 95 97.5])])




%%


fntsz=20;

c=gray(4);
c(1,:)=[];

%PSIPPS

% Create figure
figure1 = figure;

% Create axes
%axes1 = axes('Parent',figure1,...
%    'Position',[0.13 0.276080335397933 0.618426355594407 0.560273831268734]);

aboxplot3(dataplot,'colormap',c,'labels',{'CMIP6','CMIP5','low-ice EC-EARTH + 2\circC','low-ice PSIPPS + 2\circC','EC-EARTH','PSIPPS'})
ymax=round(max(max(max(dataplot))));
ylim([-50 ymax])



%RCMs+MLR uncertainties RCMs uncertainties
ylabel('Impact (percent)','FontSize',fntsz);
xlabel('Scenario','FontSize',fntsz);
legend({'Models uncertainties', 'Models+MLR uncertainties'},'Location','northoutside','Orientation','horizontal','FontSize',fntsz)
xticklabel_rotate([],45,[]);

gridxy([],[0:50:ymax],'Color','k','Linestyle',':');
set(gca,'FontSize',fntsz)

nomeout=[dir_fig,'boxplot_BA_all_' version_name ];
print( gcf, '-depsc2', nomeout ,'-painters')

