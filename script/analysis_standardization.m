clear all
close all

%%
years_model=[1981:2000 2031:2050];
years=(1951:2020);
alpha = 0.05;
NB=1000;
dir_fire='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/';
%dir_out='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/paper/figures_temp/';
dir_out='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/paper/figures_temp/';
dir_sim_pr='/Users/marco/Documents/dati/CMIP5_CMIP6/pr/ssp585/';
dir_sim_tasmax='/Users/marco/Documents/dati/CMIP5_CMIP6/tasmax/ssp585/';

domain='baileys';

%%
savename = [dir_fire,'fire_1950_2020_',domain,'_mjjas_gt300.mat'];
load(savename,'FIRE')
anni_fire=1950:2020;
[C,IA,IB] = intersect(years,anni_fire);
FIRE=FIRE(IB)/1000;

figure; plot(FIRE)

yor=(log(FIRE));
figure; plot(yor)

%% califonia nclimgrid
filename = [dir_fire,'spi4_3_nclimgrid.mat'];
load(filename);
filename = [dir_fire,'tx_4_9_nclimgrid.mat'];
load(filename);


SP=scale_1981_2000(spi4_3,years);
TX=scale_1981_2000(tx_4_9,years);

%% regression model
%[a,b]=corrcoef(SP,TX)

Xorig1 = [ones(size(SP,1),1) SP TX scale_1981_2000(years',years)  ]; %osservati


[B1,BINT,R,RINT,STATS1] = regress((yor'),Xorig1,alpha);
Y1 = (Xorig1*B1);
[a,b]=corrcoef((yor),Y1)

figure;
plot(years,(yor),'-ko','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
hold on
plot(years,(Y1),'k:o','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',6)
legend('Observed','Modelled','Location','Best')
xlabel('Years')

bootdat = [(yor)',Xorig1];
bootb1 = bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat);
bootb1 = sort(bootb1);
bootCI1 = [bootb1(25,:); bootb1(975,:)]  %Report 5% confidence intervals
B1


%% CMIP6

listGCMs_SPI = dir(fullfile(dir_sim_pr, 'pr_*_spi4.mat'));
listGCMs_TX = dir(fullfile(dir_sim_tasmax, 'tasmax*.mat'));

k=0;
for ifile=3 %1:1% length(listGCMs_SPI)
    
    load([dir_sim_pr listGCMs_SPI(ifile).name])
    load([dir_sim_tasmax listGCMs_TX(ifile).name])
    
    
    figure;plot(tx_4_9)
    
    tx_4_9_1=nanpstd(tx_4_9);
    tx_4_9_2=scale_1981_2000(tx_4_9,years_model);
    
    
    %%
    figure;plot(1981:2000, tx_4_9_1(1:20),'-o','color',[0,0,0]/255,'LineWidth',1.5)
    hold on
    plot(1981:2000, tx_4_9_2(1:20),'-s','color',[255,0,0]/255,'LineWidth',1.5)
    legend(['Old std method, change: ',num2str(nanmean(tx_4_9_1(21:40))-nanmean(tx_4_9_1(1:20)))] ,['New std method, change: ',num2str(nanmean(tx_4_9_2(21:40))-nanmean(tx_4_9_2(1:20)))],'Location','NorthWest')
    plot(2031:2050, tx_4_9_1(21:40),'-o','color',[0,0,0]/255,'LineWidth',1.5)
    plot(2031:2050, tx_4_9_2(21:40),'-s','color',[255,0,0]/255,'LineWidth',1.5)
    gridxy([],[-2:1:5],'Color','k','Linestyle',':');
    nanmean(tx_4_9_1(21:40))-nanmean(tx_4_9_1(1:20))
    nanmean(tx_4_9_2(21:40))-nanmean(tx_4_9_2(1:20))
    
    
    ylabel('TSMAX (unit less)','FontSize',18);
    xlabel('Years','FontSize',18);
    
    set(gca,'FontSize',18)
    
    nomeout=[dir_out,'cfr_std_cmip6'];
    print( gcf, '-depsc2', nomeout ,'-painters')
    
    %%
    
    nanmean(tx_4_9(21:40))-nanmean(tx_4_9(1:20))
    
    cc
    
    
    cc
    %X = [ones(size(spi4_3,1),1) scale_1981_2000(spi4_3,years_model) scale_1981_2000(tx_4_9,years_model) ];
    
    %pred_spi(:,ifile)=scale_1981_2000(spi4_3,years_model);
    %pred_tx(:,ifile)=scale_1981_2000(tx_4_9,years_model);
    pred_cmip6(:,ifile) = B1(1) + B1(2) * scale_1981_2000(spi4_3,years_model) + B1(3) * scale_1981_2000(tx_4_9,years_model);
    
    %     ctrl_ens(ifile)=nanmean(exp(pred_cmip6(1:20,ifile)));
    %     ice_ens(ifile)=nanmean(exp(pred_cmip6(21:40,ifile)));
    %
    %     changes(ifile)=100*(ice_ens(ifile)-ctrl_ens(ifile))/ctrl_ens(ifile);
    %
    %     for ib=1:1000
    %         k=k+1;
    %         %pred = bootb1(ib,3) * spi6_4 + bootb1(ib,4) * tx_6_7_2;
    %         pred_boot_cmip6(k,:) = bootb1(ib,1) + bootb1(ib,2) * scale_1981_2000(spi4_3,years) + bootb1(ib,3) * scale_1981_2000(tx_4_10,years);
    %         ctrl_ens_dum(k)=nanmean(exp(pred_dum(1:20)));
    %         ice_ens_dum(k)=nanmean(exp(pred_dum(21:40)));
    %         changes_boot(k)=100*(ice_ens_dum(k)-ctrl_ens_dum(k))/ctrl_ens_dum(k);
    %     end
end

