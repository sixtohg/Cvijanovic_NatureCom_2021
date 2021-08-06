%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all,clc,%% close all
%% cd('/home/sixto/Documents/MLToolbox_R2013/MeteoLab/'),init
workPath=[pwd '/../'];cd([workPath 'script/'])
%% workPath='C:/Users/sixto/Dropbox/ivana_ice_fires/review/';cd([workPath 'script/'])
years_model=[1981:2000 2031:2050];
years=(1951:2020);
alpha = 0.05;
NB=1000;
%% dir_fire='C:/Users/sixto/Dropbox/ivana_ice_fires/review/data/';
%% dir_out='C:/Users/sixto/Dropbox/ivana_ice_fires/review/paper/figures_temp/';
dir_fire=[workPath 'data/'];
dir_out=[workPath 'paper/figures_temp/'];
dir_sim_pr=[workPath 'data/CMIP5_CMIP6/pr/ssp585/'];
dir_sim_tasmax=[workPath 'data/CMIP5_CMIP6/tasmax/ssp585/'];
%% dir_sim_pr='C:/Users/sixto/Dropbox/Public/2020_Turco/review/CMIP5_CMIP6/pr/ssp585/';
%% dir_sim_tasmax='C:/Users/sixto/Dropbox/Public/2020_Turco/review/CMIP5_CMIP6/tasmax/ssp585/';

domain='baileys';
savename = [dir_fire,'fire_1950_2020_',domain,'_mjjas_gt300.mat'];
load(savename,'FIRE')
anni_fire=1950:2020;
[C,IA,IB] = intersect(years,anni_fire);
FIRE=FIRE(IB)/1000;

%% figure(1),subplot(1,2,1);plot(C,FIRE)
yor=(log(FIRE));
%% figure(1),subplot(1,2,2);plot(C,yor)

%% califonia nclimgrid
filename = [dir_fire,'spi4_3_nclimgrid.mat'];
load(filename);
filename = [dir_fire,'tx_4_9_nclimgrid.mat'];
load(filename);

SP=scale_1981_2000(spi4_3,years);
TX=scale_1981_2000(tx_4_9,years);

%% regression model: [a,b]=corrcoef(SP,TX)
%Xorig1 = [ones(size(SP,1),1) detrend(SP,'omitnan') detrend(TX,'omitnan') years' (years.*years)']; %osservati
Xorig1 = [ones(size(SP,1),1) SP TX years' (years.*years)']; %osservati
%Xorig1 = [ones(size(SP,1),1) SP TX years' ]; %osservati

indYears1=[1:length(years)];%% All the period
indYears2=[31:70];%% 1981:2020

[B1] = regress((yor(indYears1)'),Xorig1(indYears1,:),alpha);
[B2] = regress((yor(indYears2)'),Xorig1(indYears2,:),alpha);

Y1 = (Xorig1*B1);
Y2 = (Xorig1*B2);
%Y22 = B2(1) + B2(2) * SP + B2(3) * TX;

[a,b]=corrcoef((yor),Y1)

%% Si lo primero es calcular los cambios del modelo clima-ba con observaciones y trend, entre 2001-2020 y 1981-2020, 
%disp(sprintf('The observed (sum) difference between 2001-2020 and 1981-2000 is: %f - %f = %f', sum(exp(yor(51:end)')), sum(exp(yor(31:50)')), sum(exp(yor(51:end)'))- sum(exp(yor(31:50)'))))
%disp(sprintf('The observed (mean) difference between 2001-2020 and 1981-2000 is: %f - %f = %f', mean(exp(yor(51:end)')), mean(exp(yor(31:50)')), mean(exp(yor(51:end)')) - mean(exp(yor(31:50)'))))
disp(sprintf('The observed relative difference between 2001-2020 and 1981-2000 is: %f', 100*(mean(exp(yor(51:end)')) - mean(exp(yor(31:50)')))/mean(exp(yor(31:50)'))))
%[ci,bootstat] = bootci(NB,@(x)[mean(x) ],FIRE(31:50))
%change_obs= 100*(mean(exp(yor(51:end)')) - mean(exp(yor(31:50)')))/mean(exp(yor(31:50)'));
change_obs= 100*(mean((FIRE(51:end)')) - mean((FIRE(31:50)')))/mean((FIRE(31:50)'));
[~,bootstat1] = bootci(NB,@(x)[mean(x) ],FIRE(31:50));
[~,bootstat2] = bootci(NB,@(x)[mean(x) ],FIRE(51:70));

for ib=1:NB
    change_obs_boot(:,ib)=100*(bootstat2(ib) - bootstat1(ib))/bootstat1(ib);    
end


%% cambios del modelo con su intervalo de confianza.
bootdat1 = [(yor(indYears1))',Xorig1(indYears1,:)];
bootb1 = bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat1);

bootdat2 = [(yor(indYears2))',Xorig1(indYears2,:)];
bootb2 = bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat2);

%% bootb1 = sort(bootb1);
%% bootCI1 = [bootb1(25,:); bootb1(975,:)]; %% Report 5% confidence intervals
for ib=1:NB
    %aux_2001_2020(ib)=mean(exp(Y1_boot(51:end,ib)));
    %aux_1981_2000(ib)=mean(exp(Y1_boot(31:50,ib)));
    
    Y1_boot(:,ib) = (Xorig1*bootb1(ib,:)'); 
    Y1_boot_climate_drivers(:,ib) = bootb1(ib,1) + bootb1(ib,2) * SP + bootb1(ib,3)  * TX + bootb1(ib,4) * mean((1981:2000)')  + bootb1(ib,5) * mean((1981:2000).*(1981:2000))';   
    change_alldriver1(ib)=100*(mean(exp(Y1_boot(51:end,ib))) - mean(exp(Y1_boot(31:50,ib))))/mean(exp(Y1_boot(31:50,ib)));
    change_climate_drivers1(ib)=100*(mean(exp(Y1_boot_climate_drivers(51:end,ib))) - mean(exp(Y1_boot_climate_drivers(31:50,ib))))/mean(exp(Y1_boot_climate_drivers(31:50,ib)));
    
    Y2_boot(:,ib) = (Xorig1*bootb2(ib,:)'); 
    Y2_boot_climate_drivers(:,ib) = bootb2(ib,1) + bootb2(ib,2) * SP + bootb2(ib,3)  * TX + bootb2(ib,4) * mean((1981:2000)')  + bootb2(ib,5) * mean((1981:2000).*(1981:2000))';   
    change_alldriver2(ib)=100*(mean(exp(Y2_boot(51:end,ib))) - mean(exp(Y2_boot(31:50,ib))))/mean(exp(Y2_boot(31:50,ib)));
    change_climate_drivers2(ib)=100*(mean(exp(Y2_boot_climate_drivers(51:end,ib))) - mean(exp(Y2_boot_climate_drivers(31:50,ib))))/mean(exp(Y2_boot_climate_drivers(31:50,ib)));
    
end

%Y1_CI=prctile(Y1_boot',[2.5 97.5]);Y1_CI=Y1_CI';
%aux_2001_2020_CI=prctile(aux_2001_2020',[2.5 97.5]);aux_2001_2020_CI=aux_2001_2020_CI';
%aux_1981_2000_CI=prctile(aux_1981_2000',[2.5 97.5]);aux_1981_2000_CI=aux_1981_2000_CI';
change_alldriver_CI1=prctile(change_alldriver1',[2.5 97.5]);change_alldriver_CI1=change_alldriver_CI1';
change_climate_drivers_CI1=prctile(change_climate_drivers1',[2.5 97.5]);change_climate_drivers_CI1=change_climate_drivers_CI1';
change_alldriver_CI2=prctile(change_alldriver2',[2.5 97.5]);change_alldriver_CI2=change_alldriver_CI2';
change_climate_drivers_CI2=prctile(change_climate_drivers2',[2.5 97.5]);change_climate_drivers_CI2=change_climate_drivers_CI2';

% disp(sprintf('The modelled (sum) difference between 2001-2020 and 1981-2020 is: %f - %f = %f', sum(exp(Y1(51:end))), sum(exp(Y1(31:50))), sum(exp(Y1(51:end)))- sum(exp(Y1(31:50)))))
% disp(sprintf('The modelled (mean) difference between 2001-2020 and 1981-2020 is: %f - %f = %f', mean(exp(Y1(51:end))), mean(exp(Y1(31:50))), mean(exp(Y1(51:end)))- mean(exp(Y1(31:50)))))
% disp(sprintf('The modelled relative difference between 2001-2020 and 1981-2020 is: %f (%f, %f)', 100*(mean(exp(Y1(51:end))) - mean(exp(Y1(31:50))))/mean(exp(Y1(31:50))), 100*(mean(exp(Y1_CI(51:end,1))) - mean(exp(Y1_CI(31:50,1))))/mean(exp(Y1_CI(31:50,1))), 100*(mean(exp(Y1_CI(51:end,2))) - mean(exp(Y1_CI(31:50,2))))/mean(exp(Y1_CI(31:50,2)))))

%disp(sprintf('The modelled (mean)  2001-2020 is: %f  %f  %f', mean(exp(Y1(51:end))), aux_2001_2020_CI(1), aux_2001_2020_CI(2) ))
%disp(sprintf('The modelled (mean)  1981-2000 is: %f  %f  %f', mean(exp(Y1(31:50))), aux_1981_2000_CI(1), aux_1981_2000_CI(2) ))
disp(sprintf('The modelled (all drivers, all period) relative difference between 2001-2020 and 1981-2020 is: %f (%f, %f)', prctile(change_alldriver1',[50]), change_alldriver_CI1(1) , change_alldriver_CI1(2) ))
disp(sprintf('The modelled (all drivers, 1981-2020 period) relative difference between 2001-2020 and 1981-2020 is: %f (%f, %f)', prctile(change_alldriver2',[50]), change_alldriver_CI2(1) , change_alldriver_CI2(2) ))

disp(sprintf('The modelled (climate drivers, all period) relative difference between 2001-2020 and 1981-2020 is: %f (%f, %f)', prctile(change_climate_drivers1',[50]), change_climate_drivers_CI1(1) , change_climate_drivers_CI1(2) ))
disp(sprintf('The modelled (climate drivers, 1981-2020 period) relative difference between 2001-2020 and 1981-2020 is: %f (%f, %f)', prctile(change_climate_drivers2',[50]), change_climate_drivers_CI2(1) , change_climate_drivers_CI2(2) ))


%%

y=[change_obs prctile(change_alldriver1',[50]) prctile(change_alldriver2',[50]) prctile(change_climate_drivers1',[50]) prctile(change_climate_drivers2',[50])];
err=ones(1,5,2)*NaN;
err(1,:,1) = [change_obs-prctile(change_obs_boot',2.5) prctile(change_alldriver1',[50])-change_alldriver_CI1(1) prctile(change_alldriver2',[50])-change_alldriver_CI2(1) prctile(change_climate_drivers1',[50])-change_climate_drivers_CI1(1) prctile(change_climate_drivers2',[50])-change_climate_drivers_CI2(1)];
err(1,:,2) = [prctile(change_obs_boot',97.5)-change_obs  change_alldriver_CI1(2)-prctile(change_alldriver1',[50]) change_alldriver_CI2(2)-prctile(change_alldriver2',[50]) change_climate_drivers_CI1(2)-prctile(change_climate_drivers1',[50]) change_climate_drivers_CI2(2)-prctile(change_climate_drivers2',[50])];
% err(2,:,1) = [0 ALLD2-ALLD2low CLID2-CLID2low];
% err(2,:,2) = [0 ALLD2up-ALLD2 CLID2up-CLID2];
% err(3,:,1) = [0 ALLD3-ALLD3low CLID3-CLID3low];
% err(3,:,2) = [0 ALLD3up-ALLD3 CLID3up-CLID3];




figure(1); clf; 
%[17 17 17]
hold on
for i = 1:length(y)
    h=bar(i,y(i));
    errorbar(i, y(:,i), err(:,i,1),err(:,i,2), 'LineStyle', 'none', ... 
        'Color', [150 150 150]/255, 'LineWidth', 1);
    
    if (i == 1)
        set(h,'FaceColor','k');
    elseif (i ==2 || i == 3)
        set(h,'FaceColor','r');
    else
        set(h,'FaceColor',[255,165,0]/255);
    end
end


% % Set Axis properties
%set(gca,'xticklabel',{'1951-2000'; '1981-2000'; '20001-2020'});
%set(gca,'xticklabel',{'1951-2020'});
% ylim([200 360])
title('Burned-area changes (2001-2020 vs 1981-2000) [%]')
ylabel('Changes [%]','Fontsize',18)
%xlabel('Changes')

set(gca,'XTick',1:5);
set(gca,'XTickLabel',{'Obs','All drivers; all period','All drivers; 1981-2020 period','Climate drivers; all period','Climate drivers; 1981-2020 period'},'Fontsize',10);
htick=xticklabel_rotate([],45,[]);

%xticklabels({'Obs','All drivers','All drivers','Climate drivers','Climate drivers'})
%htick=xticklabel_rotate([],45,[]);

set(gca,'FontSize',18)
set(gcf, 'PaperPositionMode', 'auto','renderer', 'painters');
file=[dir_out,'figure1y_temp_changes.eps']
print( gcf, '-depsc2', file ,'-painters')



cc
%Y1_CI_changes = 100*(mean(exp(Y1_boot(51:end,:))) - mean(exp(Y1_boot(31:50,:))))./mean(exp(Y1_boot(31:50,:))); 
%Y1_CI_changes=prctile(Y1_CI_changes',[2.5 97.5]);
%disp(sprintf('The modelled relative difference between 2001-2020 and 1981-2020 is: %f (%f, %f)', 100*(mean(exp(Y1(51:end))) - mean(exp(Y1(31:50))))/mean(exp(Y1(31:50))), Y1_CI_changes(1), Y1_CI_changes(2)))

%% figure;
%% plot(years,(yor),'-ko','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
%% hold on
%% plot(years,(Y1),'k:o','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',6)
%% plot(years,prctile(Y1_boot',[50]),'r-')
%% legend('Observed','Modelled','Location','Best')
%% xlabel('Years')
%% plot(years,prctile(Y1_boot',[2.5]),'g--')
%% plot(years,prctile(Y1_boot',[97.5]),'g--')

%% figure;
%% plot([1 2 3],[mean(exp(yor(1:50)')) mean(exp(yor(31:50)')) mean(exp(yor(51:end)'))],'k.','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
%% hold on, plot([1 2 3],[mean(exp(Y1(1:50)')) mean(exp(Y1(31:50)')) mean(exp(Y1(51:end)'))],'r.','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6)

%% mira, lo que m??s me preocupas es que no tenemos modelos, creo, que van bien para responder a este revisor:
%% 
%% " I would suggest the author add one panel in figure 1 to show the decadal burned-area changes (2000-2019 vs. 1980-1999) 
%% between observed data and the regression wildfire model."
%% 
%% puedes hacer unas pruebas con distintos periodos (1951-2020 o 1981-2020) y modelos (t, o t^2) para encontrar un modelo que 
%% vaya bien para decadal burned-area changes (2001-2021 vs. 1981-2000)?
%% 
%% y si no hay, uno que vaya bien para los valores de ba en los dos periodos, 2001-2021 y 1981-2000, que igual los "changes" como son % se amplifican los errores.
auxBarplot=repmat(NaN,4,7);
for i=1:7
	%% auxBarplot(:,i)=[mean(exp(yor((i-1)*10+1:i*10)')) mean(exp(Y1((i-1)*10+1:i*10)')) prctile(mean(exp(Y1_CI((i-1)*10+1:i*10,1))) mean(exp(Y1_CI((i-1)*10+1:i*10,2)))];
	auxBarplot(:,i)=[mean(exp(yor((i-1)*10+1:i*10)')) mean(exp(Y1((i-1)*10+1:i*10)')) prctile(mean(exp(Y1_boot((i-1)*10+1:i*10,:)))',2.5) prctile(mean(exp(Y1_boot((i-1)*10+1:i*10,:)))',97.5)];
end
figure,bar([1:7],auxBarplot,'grouped')

figure,bar([1:6],[mean(exp(yor(1:50)')) mean(exp(Y1(1:50)')) mean(exp(yor(31:50)')) mean(exp(Y1(31:50)')) mean(exp(yor(51:end)')) mean(exp(Y1(51:end)'))])
hold on, plot([2 4 6],[mean(exp(Y1_CI(1:50,1))) mean(exp(Y1_CI(31:50,1))) mean(exp(Y1_CI(51:end,1)))],'g.','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6)
hold on, plot([2 4 6],[mean(exp(Y1_CI(1:50,2))) mean(exp(Y1_CI(31:50,2))) mean(exp(Y1_CI(51:end,2)))],'g.','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6)

%% CMIP6
listGCMs_SPI = dir(fullfile(dir_sim_pr, 'pr_*_spi4.mat'));
listGCMs_TX = dir(fullfile(dir_sim_tasmax, 'tasmax*.mat'));
noSPI=0;
noTX=0;
k=0;
changes_cmip6=repmat(NaN,length(listGCMs_SPI),6);
clim_cmip6=repmat(NaN,length(listGCMs_SPI),6);
for ifile=1:length(listGCMs_SPI)
    if isempty(strfind(listGCMs_SPI(ifile).name,'CMIP6_MIROC6'))
		load([dir_sim_pr listGCMs_SPI(ifile).name]),if noSPI==1,spi4_3(21:40)=nanmean(spi4_3(1:20));end
		load([dir_sim_tasmax listGCMs_TX(ifile).name]),if noTX==1,tx_4_9(21:40)=nanmean(tx_4_9(1:20));end
		pred_cmip6(:,ifile) = B1(1) + B1(2) * scale_1981_2000(spi4_3,years_model) + B1(3) * scale_1981_2000(tx_4_9,years_model);
		ctrl_ens=nanmean(exp(pred_cmip6(1:20,ifile)-B1(1)));%% Due to the first coefficient the model explodes and it is simplified in the relative delta change.
		ice_ens=nanmean(exp(pred_cmip6(21:40,ifile)-B1(1)));
		changes_cmip6(ifile,1)=100*(ice_ens-ctrl_ens)/ctrl_ens;
		changes_cmip6(ifile,2)=(ice_ens-ctrl_ens)*exp(B1(1));
		changes_cmip6(ifile,3)=nanmean(tx_4_9(21:40))-nanmean(tx_4_9(1:20));
		changes_cmip6(ifile,4)=nanmean(spi4_3(21:40))-nanmean(spi4_3(1:20));
		Y1_boot=repmat(NaN,length(years_model),NB);
		for ib=1:NB
			Y1_boot(:,ib) = exp(bootb1(ib,2) * scale_1981_2000(spi4_3,years_model) + bootb1(ib,3) * scale_1981_2000(tx_4_9,years_model));%% Due to the first coefficient the model explodes and it is simplified in the relative delta change.
		end
		aux_ens=(nanmean(Y1_boot(21:40,:)) - nanmean(Y1_boot(1:20,:)))./nanmean(Y1_boot(1:20,:));
		Y1_CI=100*prctile(aux_ens',[2.5 97.5]);
		changes_cmip6(ifile,5:6)=Y1_CI;
		Y1_boot=repmat(NaN,length(years_model),NB);
		for ib=1:NB
			Y1_boot(:,ib) = exp(bootb1(ib,1) + bootb1(ib,2) * scale_1981_2000(spi4_3,years_model) + bootb1(ib,3) * scale_1981_2000(tx_4_9,years_model) + bootb1(ib,4) *years_model(:) + bootb1(ib,5) *years_model(:).^2);
		end
		aux_ens=nanmean(Y1_boot(1:20,:));
		Y1_CI=prctile(aux_ens',[2.5 97.5]);
		aux_cmip6 = exp(B1(1) + B1(2) * scale_1981_2000(spi4_3,years_model) + B1(3) * scale_1981_2000(tx_4_9,years_model) + B1(4)*years_model(:) + B1(5)*years_model(:).^2);
		clim_cmip6(ifile,1:3)=[nanmean(aux_cmip6(1:20)) Y1_CI(:)'];
	end
end
clim_cmip6(:,4:6)=[clim_cmip6(:,1).*(changes_cmip6(:,1)+100)/100 clim_cmip6(:,1).*(changes_cmip6(:,5)+100)/100 clim_cmip6(:,1).*(changes_cmip6(:,6)+100)/100];

%% CMIP5
dir_sim_pr=[workPath 'data/CMIP5_CMIP6/pr/rcp85/'];
dir_sim_tasmax=[workPath 'data/CMIP5_CMIP6/tasmax/rcp85/'];
listGCMs={'ACCESS1-3','ACCESS1-3','BCC-CSM1-1','CNRM-CM5','CNRM-CM5','CanESM2','EC-EARTH','EC-EARTH','FGOALS-g2','GFDL-ESM2G', ...
    'GFDL-ESM2M','inmcm4','inmcm4','IPSL-CM5A-LR','MIROC-ESM','MIROC5','MPI-ESM-LR','MRI-CGCM3','NorESM1-M','NorESM1-M','HadGEM2-ES'};
k=0;
changes_cmip5=repmat(NaN,length(listGCMs),6);
clim_cmip5=repmat(NaN,length(listGCMs),6);
refTmax=repmat(NaN,20,length(listGCMs_SPI));
refSPI=repmat(NaN,20,length(listGCMs_SPI));
for ifile=1:length(listGCMs)
	load([dir_sim_pr,'spi4_3_',num2str(ifile-1),'_nostd_',domain,'.mat'])
	load([dir_sim_tasmax,'tx_4_9_',num2str(ifile-1),'_nostd_',domain,'.mat'])
	if noSPI==1,spi4_3(21:40)=nanmean(spi4_3(1:20));end
	if noTX==1,tx_4_9(21:40)=nanmean(tx_4_9(1:20));end
	pred_cmip5(:,ifile) = B1(1) + B1(2) * scale_1981_2000(spi4_3,years_model) + B1(3) * scale_1981_2000(tx_4_9,years_model);
	ctrl_ens=nanmean(exp(pred_cmip5(1:20,ifile)-B1(1)));
	ice_ens=nanmean(exp(pred_cmip5(21:40,ifile)-B1(1)));
	changes_cmip5(ifile,1)=100*(ice_ens-ctrl_ens)/ctrl_ens;
	changes_cmip5(ifile,2)=(ice_ens-ctrl_ens)*exp(B1(1));
	changes_cmip5(ifile,3)=nanmean(tx_4_9(21:40))-nanmean(tx_4_9(1:20));
	changes_cmip5(ifile,4)=nanmean(spi4_3(21:40))-nanmean(spi4_3(1:20));
	refTmax(:,ifile)=tx_4_9(1:20);refSPI(:,ifile)=spi4_3(1:20);
	Y1_boot=repmat(NaN,length(years_model),NB);
	for ib=1:NB
		Y1_boot(:,ib) = exp(bootb1(ib,2) * scale_1981_2000(spi4_3,years_model) + bootb1(ib,3) * scale_1981_2000(tx_4_9,years_model));
	end
	aux_ens=(nanmean(Y1_boot(21:40,:)) - nanmean(Y1_boot(1:20,:)))./nanmean(Y1_boot(1:20,:));
	Y1_CI=100*prctile(aux_ens',[2.5 97.5]);
	changes_cmip5(ifile,5:6)=Y1_CI;
	Y1_boot=repmat(NaN,length(years_model),NB);
	for ib=1:NB
		Y1_boot(:,ib) = exp(bootb1(ib,1) + bootb1(ib,2) * scale_1981_2000(spi4_3,years_model) + bootb1(ib,3) * scale_1981_2000(tx_4_9,years_model) + bootb1(ib,4) *years_model(:) + bootb1(ib,5) *years_model(:).^2);
	end
	aux_ens=nanmean(Y1_boot(1:20,:));
	Y1_CI=prctile(aux_ens',[2.5 97.5]);
	aux_cmip5 = exp(B1(1) + B1(2) * scale_1981_2000(spi4_3,years_model) + B1(3) * scale_1981_2000(tx_4_9,years_model) + B1(4)*years_model(:) + B1(5)*years_model(:).^2);
	clim_cmip5(ifile,1:3)=[nanmean(aux_cmip5(1:20)) Y1_CI(:)'];
end
clim_cmip5(:,4:6)=[clim_cmip5(:,1).*(changes_cmip5(:,1)+100)/100 clim_cmip5(:,1).*(changes_cmip5(:,5)+100)/100 clim_cmip5(:,1).*(changes_cmip5(:,6)+100)/100];

refTmax=nanmean(refTmax,2);refSPI=nanmean(refSPI,2);
auxP=exp(B1(1) + B1(2)*refSPI+ B1(3)*refTmax);
changes_spi=-1.0:0.1:1.0;changes_tmax=-2.0:0.20:4.0;
[a1,a2]=meshgrid(changes_spi,changes_tmax);
changes_neutral=repmat(NaN,size(a1));
for i=1:size(a1,1)
	for j=1:size(a1,2)
		auxF=exp(B1(2) * (refSPI + a1(i,j)) + B1(3) * (refTmax + a2(i,j)));
		changes_neutral(i,j)=100*(nanmean(auxF)-nanmean(auxP))/nanmean(auxP);
	end
end
changes_neutral=changes_neutral';

%% low-ice
%% dir_sim='C:/Users/sixto/Dropbox/ivana_ice_fires/review/data/low_ice/';
dir_sim=[workPath 'data/low_ice/'];
k=0;changes_lowice=repmat(NaN,12,6);
clim_lowice=repmat(NaN,12,6);
for ifile=1:12
    load([dir_sim 'spi4_3_' num2str(ifile-1) '_nostd_' char(domain) '.mat'])
    load([dir_sim  'tx_4_9_' num2str(ifile-1) '_nostd_' char(domain) '.mat'])
	if noSPI==1,spi4_3(21:40)=nanmean(spi4_3(1:20));end
	if noTX==1,tx_4_9(21:40)=nanmean(tx_4_9(1:20));end
    pred_lowice(:,ifile) = B1(1) + B1(2) * scale_1981_2000(spi4_3,years_model) + B1(3) * scale_1981_2000(tx_4_9,years_model);
	ctrl_ens=nanmean(exp(pred_lowice(1:20,ifile)-B1(1)));
	ice_ens=nanmean(exp(pred_lowice(21:40,ifile)-B1(1)));
	changes_lowice(ifile,1)=100*(ice_ens-ctrl_ens)/ctrl_ens;
	changes_lowice(ifile,2)=(ice_ens-ctrl_ens)*exp(B1(1));
	changes_lowice(ifile,3)=nanmean(tx_4_9(21:40))-nanmean(tx_4_9(1:20));
	changes_lowice(ifile,4)=nanmean(spi4_3(21:40))-nanmean(spi4_3(1:20));
	
	Y1_boot=repmat(NaN,length(years_model),NB);
	for ib=1:NB
		Y1_boot(:,ib) = exp(bootb1(ib,2) * scale_1981_2000(spi4_3,years_model) + bootb1(ib,3) * scale_1981_2000(tx_4_9,years_model));
	end
	aux_ens=(nanmean(Y1_boot(21:40,:)) - nanmean(Y1_boot(1:20,:)))./nanmean(Y1_boot(1:20,:));
	Y1_CI=100*prctile(aux_ens',[2.5 97.5]);
	changes_lowice(ifile,5:6)=Y1_CI;
	Y1_boot=repmat(NaN,length(years_model),NB);
	for ib=1:NB
		Y1_boot(:,ib) = exp(bootb1(ib,1) + bootb1(ib,2) * scale_1981_2000(spi4_3,years_model) + bootb1(ib,3) * scale_1981_2000(tx_4_9,years_model) + bootb1(ib,4) *years_model(:) + bootb1(ib,5) *years_model(:).^2);
	end
	aux_ens=nanmean(Y1_boot(1:20,:));
	Y1_CI=prctile(aux_ens',[2.5 97.5]);
	aux_cmip5 = exp(B1(1) + B1(2) * scale_1981_2000(spi4_3,years_model) + B1(3) * scale_1981_2000(tx_4_9,years_model) + B1(4)*years_model(:) + B1(5)*years_model(:).^2);
	clim_lowice(ifile,1:3)=[nanmean(aux_cmip5(1:20)) Y1_CI(:)'];
end
clim_lowice(:,4:6)=[clim_lowice(:,1).*(changes_lowice(:,1)+100)/100 clim_lowice(:,1).*(changes_lowice(:,5)+100)/100 clim_lowice(:,1).*(changes_lowice(:,6)+100)/100];

%% low ice ec-earth
%% dir_sim='C:/Users/sixto/Dropbox/ivana_ice_fires/review/data/low_ice_ecearth/';
dir_sim=[workPath 'data/low_ice_ecearth/'];
k=0;changes_EC=repmat(NaN,10,6);
clim_EC=repmat(NaN,10,6);
for ifile=1:10
    load([dir_sim 'spi4_3_' num2str(ifile-1) '_nostd_' char(domain) '.mat'])
    load([dir_sim  'tx_4_9_' num2str(ifile-1) '_nostd_' char(domain) '.mat'])
	if noSPI==1,spi4_3(21:40)=nanmean(spi4_3(1:20));end
	if noTX==1,tx_4_9(21:40)=nanmean(tx_4_9(1:20));end
    pred_ecearth(:,ifile) = B1(1) + B1(2) * scale_1981_2000(spi4_3,years_model) + B1(3) * scale_1981_2000(tx_4_9,years_model);
	ctrl_ens=nanmean(exp(pred_ecearth(1:20,ifile)-B1(1)));
	ice_ens=nanmean(exp(pred_ecearth(21:40,ifile)-B1(1)));
	changes_EC(ifile,1)=100*(ice_ens-ctrl_ens)/ctrl_ens;
	changes_EC(ifile,2)=(ice_ens-ctrl_ens)*exp(B1(1));
	changes_EC(ifile,3)=nanmean(tx_4_9(21:40))-nanmean(tx_4_9(1:20));
	changes_EC(ifile,4)=nanmean(spi4_3(21:40))-nanmean(spi4_3(1:20));
	
	Y1_boot=repmat(NaN,length(years_model),NB);
	for ib=1:NB
		Y1_boot(:,ib) = exp(bootb1(ib,2) * scale_1981_2000(spi4_3,years_model) + bootb1(ib,3) * scale_1981_2000(tx_4_9,years_model));
	end
	aux_ens=(nanmean(Y1_boot(21:40,:)) - nanmean(Y1_boot(1:20,:)))./nanmean(Y1_boot(1:20,:));
	Y1_CI=100*prctile(aux_ens',[2.5 97.5]);
	changes_EC(ifile,5:6)=Y1_CI;
	Y1_boot=repmat(NaN,length(years_model),NB);
	for ib=1:NB
		Y1_boot(:,ib) = exp(bootb1(ib,1) + bootb1(ib,2) * scale_1981_2000(spi4_3,years_model) + bootb1(ib,3) * scale_1981_2000(tx_4_9,years_model) + bootb1(ib,4) *years_model(:) + bootb1(ib,5) *years_model(:).^2);
	end
	aux_ens=nanmean(Y1_boot(1:20,:));
	Y1_CI=prctile(aux_ens',[2.5 97.5]);
	aux_cmip5 = exp(B1(1) + B1(2) * scale_1981_2000(spi4_3,years_model) + B1(3) * scale_1981_2000(tx_4_9,years_model) + B1(4)*years_model(:) + B1(5)*years_model(:).^2);
	clim_EC(ifile,1:3)=[nanmean(aux_cmip5(1:20)) Y1_CI(:)'];
end
clim_EC(:,4:6)=[clim_EC(:,1).*(changes_EC(:,1)+100)/100 clim_EC(:,1).*(changes_EC(:,5)+100)/100 clim_EC(:,1).*(changes_EC(:,6)+100)/100];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Estimation of the contribution of each component %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sep_cmip6=[100*(exp(B1(2)*changes_cmip6(:,4) + B1(3)*changes_cmip6(:,3)) - 1) 100*(exp(B1(2)*changes_cmip6(:,4) + 0*changes_cmip6(:,3)) - 1) 100*(exp(0*changes_cmip6(:,4) + B1(3)*changes_cmip6(:,3)) - 1)];
sep_cmip5=[100*(exp(B1(2)*changes_cmip5(:,4) + B1(3)*changes_cmip5(:,3)) - 1) 100*(exp(B1(2)*changes_cmip5(:,4) + 0*changes_cmip5(:,3)) - 1) 100*(exp(0*changes_cmip5(:,4) + B1(3)*changes_cmip5(:,3)) - 1)];
sep_lowice=[100*(exp(B1(2)*changes_lowice(:,4) + B1(3)*changes_lowice(:,3)) - 1) 100*(exp(B1(2)*changes_lowice(:,4) + 0*changes_lowice(:,3)) - 1) 100*(exp(0*changes_lowice(:,4) + B1(3)*changes_lowice(:,3)) - 1)];
sep_EC=[100*(exp(B1(2)*changes_EC(:,4) + B1(3)*changes_EC(:,3)) - 1) 100*(exp(B1(2)*changes_EC(:,4) + 0*changes_EC(:,3)) - 1) 100*(exp(0*changes_EC(:,4) + B1(3)*changes_EC(:,3)) - 1)];

disp('Low-Ice:')
disp([min(sep_lowice(:,1)) max(sep_lowice(:,1)) prctile(sep_lowice(:,1),[2.5 5 25 50 75 95 97.5])]);
disp([min(sep_lowice(:,2)) max(sep_lowice(:,2)) prctile(sep_lowice(:,2),[2.5 5 25 50 75 95 97.5])]);
disp([min(sep_lowice(:,3)) max(sep_lowice(:,3)) prctile(sep_lowice(:,3),[2.5 5 25 50 75 95 97.5])]);
disp('Low-Ice EC-Earth:')
disp([min(sep_EC(:,1)) max(sep_EC(:,1)) prctile(sep_EC(:,1),[2.5 5 25 50 75 95 97.5])]);
disp([min(sep_EC(:,2)) max(sep_EC(:,2)) prctile(sep_EC(:,2),[2.5 5 25 50 75 95 97.5])]);
disp([min(sep_EC(:,3)) max(sep_EC(:,3)) prctile(sep_EC(:,3),[2.5 5 25 50 75 95 97.5])]);
disp('CMIP5:')
disp([min(sep_cmip5(:,1)) max(sep_cmip5(:,1)) prctile(sep_cmip5(:,1),[2.5 5 25 50 75 95 97.5])]);
disp([min(sep_cmip5(:,2)) max(sep_cmip5(:,2)) prctile(sep_cmip5(:,2),[2.5 5 25 50 75 95 97.5])]);
disp([min(sep_cmip5(:,3)) max(sep_cmip5(:,3)) prctile(sep_cmip5(:,3),[2.5 5 25 50 75 95 97.5])]);
disp('CMIP6:')
disp([min(sep_cmip6(:,1)) max(sep_cmip6(:,1)) prctile(sep_cmip6(:,1),[2.5 5 25 50 75 95 97.5])]);
disp([min(sep_cmip6(:,2)) max(sep_cmip6(:,2)) prctile(sep_cmip6(:,2),[2.5 5 25 50 75 95 97.5])]);
disp([min(sep_cmip6(:,3)) max(sep_cmip6(:,3)) prctile(sep_cmip6(:,3),[2.5 5 25 50 75 95 97.5])]);

%% Low-Ice (min, max, percentiles 2.5 5 25 50 75 95 97.5): Only including time
%% All  -35.0994   80.2047  -35.0994  -33.4083  -10.0109    5.3273   36.6668   77.9028   80.2047
%% SPI   -2.8155   13.3796   -2.8155   -2.7782   -0.6340    4.4697    6.7785   13.3475   13.3796
%% TX   -42.7581   84.7167  -42.7581  -40.8134  -11.6415    1.6302   31.7666   81.8361   84.7167
%% 
%% Low-Ice EC-Earth:
%% All   -8.2828   37.2927   -8.2828   -8.2828    1.3841    9.6377   31.3828   37.2927   37.2927
%% SPI   -8.2283   10.0623   -8.2283   -8.2283   -0.7602    2.5795    9.2692   10.0623   10.0623
%% TX    -9.0869   32.9055   -9.0869   -9.0869    0.5015   10.9745   24.7409   32.9055   32.9055
%% 
%% CMIP5:
%% All  122.5789  568.9350  122.7201  125.6844  266.9017  334.4812  468.9011  562.2183  568.6297
%% SPI  -11.9061   12.0545  -11.8972  -11.7096   -7.0994   -4.0065   -0.7606   11.7972   12.0428
%% TX   144.4998  581.4379  144.4998  144.4998  263.2018  339.4233  484.9319  568.8470  580.8656
%% 
%% CMIP6:
%% All     195.3    1286.2     195.3     195.8     325.2     475.2     619.3    1139.2    1273.9
%% SPI  -15.0554   15.1554  -14.8948  -13.1280   -6.5459   -2.7141    2.3832    8.7866   14.6247
%% TX      185.2    1330.2     188.4     223.3     365.9     470.4     636.6    1248.9    1323.4

%% Low-Ice (min, max, percentiles 2.5 5 25 50 75 95 97.5): Including time and time^2:
%% All  -15.1508   52.3381  -15.1508  -14.8283   -6.7136   10.7949   38.7667   51.6294   52.3381
%% SPI   -5.8705   30.4747   -5.8705   -5.7939   -1.3133    9.7168   14.9054   30.3968   30.4747
%% TX   -34.9689   60.5316  -34.9689  -33.3233   -9.1167    1.2333   23.6944   58.5638   60.5316
%% 
%% Low-Ice EC-Earth:
%% All  -12.2931   45.2981  -12.2931  -12.2931   -0.3201   13.6720   32.0133   45.2981   45.2981
%% SPI  -16.6314   22.5200  -16.6314  -16.6314   -1.6036    5.5509   20.6573   22.5200   22.5200
%% TX    -7.0846   24.5348   -7.0846   -7.0846    0.3866    8.3609   18.5914   24.5348   24.5348
%% 
%% CMIP5:
%% All   63.3302  334.8343   63.5528   68.2262  158.5528  190.8240  233.8507  329.6676  334.5995
%% SPI  -23.5508   27.2655  -23.5343  -23.1887  -14.4370   -8.2975   -1.6034   26.6480   27.2375
%% TX    99.2900  339.3749   99.2900   99.2900  170.3893  213.2320  290.5477  333.0887  339.0891
%% 
%% CMIP6:
%% All   85.4297  628.4599   88.2315  119.0511  188.2785  284.0562  333.1671  579.3258  624.3654
%% SPI  -29.2248   34.8420  -28.9353  -25.7508  -13.3606   -5.6607    5.1158   19.8550   33.5930
%% TX   124.4458  678.3532  126.3234  146.9771  227.6860  283.0267  366.5737  643.8195  675.4754
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Low-Ice:')
disp([min(changes_lowice(:,1)) max(changes_lowice(:,1)) prctile(changes_lowice(:,1),[2.5 5 25 50 75 95 97.5])]);
%% disp([min(changes_lowice(:,2)) max(changes_lowice(:,2)) prctile(changes_lowice(:,2),[5 25 50 75 95])]);
disp([min(changes_lowice(:,3)) max(changes_lowice(:,3)) prctile(changes_lowice(:,3),[2.5 5 25 50 75 95 97.5])]);
disp([min(changes_lowice(:,4)) max(changes_lowice(:,4)) prctile(changes_lowice(:,4),[2.5 5 25 50 75 95 97.5])]);
disp('Low-Ice EC-Earth:')
disp([min(changes_EC(:,1)) max(changes_EC(:,1)) prctile(changes_EC(:,1),[2.5 5 25 50 75 95 97.5])]);
%% disp([min(changes_EC(:,2)) max(changes_EC(:,2)) prctile(changes_EC(:,2),[5 25 50 75 95])]);
disp([min(changes_EC(:,3)) max(changes_EC(:,3)) prctile(changes_EC(:,3),[2.5 5 25 50 75 95 97.5])]);
disp([min(changes_EC(:,4)) max(changes_EC(:,4)) prctile(changes_EC(:,4),[2.5 5 25 50 75 95 97.5])]);
disp('CMIP5:')
disp([min(changes_cmip5(:,1)) max(changes_cmip5(:,1)) prctile(changes_cmip5(:,1),[2.5 5 25 50 75 95 97.5])]);
%% disp([min(changes_cmip5(:,2)) max(changes_cmip5(:,2)) prctile(changes_cmip5(:,2),[5 25 50 75 95])]);
disp([min(changes_cmip5(:,3)) max(changes_cmip5(:,3)) prctile(changes_cmip5(:,3),[2.5 5 25 50 75 95 97.5])]);
disp([min(changes_cmip5(:,4)) max(changes_cmip5(:,4)) prctile(changes_cmip5(:,4),[2.5 5 25 50 75 95 97.5])]);
disp('CMIP6:')
disp([min(changes_cmip6(:,1)) max(changes_cmip6(:,1)) prctile(changes_cmip6(:,1),[2.5 5 25 50 75 95 97.5])]);
%% disp([min(changes_cmip6(:,2)) max(changes_cmip6(:,2)) prctile(changes_cmip6(:,2),[5 25 50 75 95])]);
disp([min(changes_cmip6(:,3)) max(changes_cmip6(:,3)) prctile(changes_cmip6(:,3),[2.5 5 25 50 75 95 97.5])]);
disp([min(changes_cmip6(:,4)) max(changes_cmip6(:,4)) prctile(changes_cmip6(:,4),[2.5 5 25 50 75 95 97.5])]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Low-Ice (min, max, percentiles 2.5 5 25 50 75 95 97.5): Only including time
%% BA   -25.9836  105.2413  -25.9836  -24.5039    6.0190   30.3270   69.4518  103.4338  105.2413
%% TX    -0.7608    0.8368   -0.7608   -0.7209   -0.1698    0.0204    0.3751    0.8137    0.8368
%% SPI   -0.6075    0.1382   -0.6075   -0.6062   -0.3173   -0.2113    0.0313    0.1363    0.1382
%% 
%% Low-Ice EC-Earth:
%% BA   -13.8937   74.8365  -13.8937  -13.8937   -4.1665   17.5965   34.4517   74.8365   74.8365
%% TX    -0.1299    0.3879   -0.1299   -0.1299    0.0068    0.1418    0.3015    0.3879    0.3879
%% SPI   -0.4639    0.4154   -0.4639   -0.4639   -0.4289   -0.1231    0.0369    0.4154    0.4154
%% 
%% CMIP5:
%% BA    86.7230  776.6463   86.8976   90.5640  214.7810  361.2872  429.5252  770.5257  776.3681
%% TX     1.2192    2.6169    1.2192    1.2192    1.7577    2.0186    2.4085    2.5912    2.6157
%% SPI   -0.5506    0.6133   -0.5501   -0.5395    0.0370    0.1978    0.3564    0.6025    0.6128
%% 
%% CMIP6:
%% BA       90.5    1192.5      91.5     102.8     239.8     429.5     801.6    1065.6    1181.9
%% TX     1.4293    3.6278    1.4430    1.5934    2.0983    2.3740    2.7231    3.5463    3.6210
%% SPI   -0.6827    0.7894   -0.6593   -0.4020   -0.1139    0.1332    0.3275    0.6817    0.7805
%% 
%% Low-Ice (min, max, percentiles 2.5 5 25 50 75 95 97.5): Including time and time^2:
%% BA   -16.8733   90.8590  -16.8733  -16.3032    5.8068   27.0192   64.1613   89.3441   90.8590
%% TX    -0.7608    0.8368   -0.7608   -0.7209   -0.1698    0.0204    0.3751    0.8137    0.8368
%% SPI   -0.6075    0.1382   -0.6075   -0.6062   -0.3173   -0.2113    0.0313    0.1363    0.1382
%% 
%% Low-Ice EC-Earth:
%% BA   -10.1820   71.1147  -10.1820  -10.1820   -4.8139   16.1052   29.9571   71.1147   71.1147
%% TX    -0.1299    0.3879   -0.1299   -0.1299    0.0068    0.1418    0.3015    0.3879    0.3879
%% SPI   -0.4639    0.4154   -0.4639   -0.4639   -0.4289   -0.1231    0.0369    0.4154    0.4154
%% 
%% CMIP5:
%% BA    64.6039  542.0518   64.7954   68.8178  165.9628  273.9861  317.4076  535.9027  541.7723
%% TX     1.2192    2.6169    1.2192    1.2192    1.7577    2.0186    2.4085    2.5912    2.6157
%% SPI   -0.5506    0.6133   -0.5501   -0.5395    0.0370    0.1978    0.3564    0.6025    0.6128
%% 
%% CMIP6:
%% BA    62.8704  795.5009   64.0564   77.1025  187.9033  310.0049  555.1425  696.5446  787.2545
%% TX     1.4293    3.6278    1.4430    1.5934    2.0983    2.3740    2.7231    3.5463    3.6210
%% SPI   -0.6827    0.7894   -0.6593   -0.4020   -0.1139    0.1332    0.3275    0.6817    0.7805
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataPlot=repmat(NaN,23,4);
dataPlot(:,1)=changes_cmip6(:,1);dataPlot(1:size(changes_cmip5,1),2)=changes_cmip5(:,1);
dataPlot(1:size(changes_lowice,1),3)=changes_lowice(:,1);dataPlot(1:size(changes_EC,1),4)=changes_EC(:,1);
figure,boxplot(dataPlot)
set(gca, 'XTick', [1:4])
set(gca,'XTickLabel',{'CMIP6','CMIP5','Low-Ice','EC-EARTH-Low-Ice'})
ylabel('Burned Area (%)','Fontsize',10,'FontWeight','bold');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numcolors=21;
load([workPath 'data/CustomColormap.mat']);
fall=CustomColormap;clear CustomColormap
%% fall=flipud(cbrewer('div', 'PuOr', 41));
%% fall=[fall(11,:);fall(21,:);fall(23:end,:)];%% fall(1,:)=[0.1765 0 0.2941];
auxValues=linspace(-40,760,numcolors);auxAxis=(auxValues(end)-auxValues(1))/(numcolors-1);
fall(find(fall < 0))=0;fall(find(fall > 1))=1;
figure;h=surface(changes_tmax,changes_spi,changes_neutral);
axis square
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
caxis([-60 780]);colormap(fall);h=colorbar;hold on
set(h, 'YTick', [-40:40:760])
set(h,'YTickLabel',[-40:40:760])
title(h,'Changes in BA','FontWeight','bold')

set(gca,'Fontsize',12);
set(findall(gcf,'-property','FontSize'),'FontSize',12)
ylabel('Changes in SPI','Fontsize',10,'FontWeight','bold');
xlabel('Changes in TSMAX (\circC)','Fontsize',10,'FontWeight','bold');

figure,
hold on,scatter(changes_lowice(:,3),changes_lowice(:,4),100,'Marker','^','MarkerEdgeColor','black')
hold on,scatter(changes_EC(:,3),changes_EC(:,4),100,'Marker','v','MarkerEdgeColor','black')
scatter(changes_cmip5(:,3),changes_cmip5(:,4),100,'Marker','o','MarkerEdgeColor','black')
hold on,scatter(changes_cmip6(:,3),changes_cmip6(:,4),100,'Marker','s','MarkerEdgeColor','black')
caxis([-60 780]);colormap(fall);h=colorbar;hold on
legend({'Low-Ice','Low-Ice EC-Earth','CMIP5','CMIP6'},'Location','Best')

set(gca,'Fontsize',12,'xlim',[-2 4],'ylim',[-1.0 1.0]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)
ylabel('Changes in SPI','Fontsize',10,'FontWeight','bold');
xlabel('Changes in TSMAX (\circC)','Fontsize',10,'FontWeight','bold');

figure,
hold on,scatter(changes_lowice(:,3),changes_lowice(:,4),100,changes_lowice(:,1),'filled','Marker','^','MarkerEdgeColor','black')
hold on,scatter(changes_EC(:,3),changes_EC(:,4),100,changes_EC(:,1),'filled','Marker','v','MarkerEdgeColor','black')
scatter(changes_cmip5(:,3),changes_cmip5(:,4),100,changes_cmip5(:,1),'filled','Marker','o','MarkerEdgeColor','black')
hold on,scatter(changes_cmip6(:,3),changes_cmip6(:,4),100,changes_cmip6(:,1),'filled','Marker','s','MarkerEdgeColor','black')
caxis([-60 780]);colormap(fall);h=colorbar;hold on
legend({['Low-Ice (' num2str(median(changes_lowice(:,1)),'%5.1f') ')'],['Low-Ice EC-Earth (' num2str(median(changes_EC(:,1)),'%5.1f') ')'],['CMIP5 (' num2str(median(changes_cmip5(:,1)),'%5.1f') ')'],['CMIP6 (' num2str(median(changes_cmip6(:,1)),'%5.1f') ')']},'Location','Best')
%% grid on
set(h, 'YTick', [-40:40:760])
set(h,'YTickLabel',[-40:40:760])
title(h,'Changes in BA','FontWeight','bold')
set(gca,'Fontsize',12,'xlim',[-2 4],'ylim',[-1.0 1.0]);
set(findall(gcf,'-property','FontSize'),'FontSize',12)
ylabel('Changes in SPI','Fontsize',10,'FontWeight','bold');
xlabel('Changes in TSMAX (\circC)','Fontsize',10,'FontWeight','bold');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For the calculations, I think it is still missing the values in ba/year 
%% for the observation and model with observation for the two 20 years periods, 
%% and for the gcms for the 1981,2000 and 2031,2050 periods.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pred_cmip6(:,17)=NaN;pred_cmip6=pred_cmip6+B1(4)*repmat(years_model(:),1,size(pred_cmip6,2));
fid=fopen([workPath 'data/BA_cmip6.csv'],'w');
fprintf(fid,'YYYY');
for i=1:size(pred_cmip6,2)
	fprintf(fid,', GCM%02d-CMIP6',i);
end
fprintf(fid,'\r\n');
for n=1:size(pred_cmip6,1)
	fprintf(fid,'%04d',years_model(n));
	for i=1:size(pred_cmip6,2)
		fprintf(fid,', %f',exp(pred_cmip6(n,i)));
	end
	fprintf(fid,'\r\n');
end
fprintf(fid,'%04d-%04d',years_model(20),years_model(1));
for i=1:size(pred_cmip6,2)
	fprintf(fid,', %f',nanmean(exp(pred_cmip6(1:20,i))));
end
fprintf(fid,'\r\n');
fprintf(fid,'%04d-%04d',years_model(40),years_model(21));
for i=1:size(pred_cmip6,2)
	fprintf(fid,', %f',nanmean(exp(pred_cmip6(21:40,i))));
end
fprintf(fid,'\r\n');
fclose(fid);
%%%%%%%%%%%%
pred_cmip5=pred_cmip5+B1(4)*repmat(years_model(:),1,size(pred_cmip5,2));
fid=fopen([workPath 'data/BA_cmip5.csv'],'w');
fprintf(fid,'YYYY');
for i=1:size(pred_cmip5,2)
	fprintf(fid,', GCM%02d-CMIP5',i);
end
fprintf(fid,'\r\n');
for n=1:size(pred_cmip5,1)
	fprintf(fid,'%04d',years_model(n));
	for i=1:size(pred_cmip5,2)
		fprintf(fid,', %f',exp(pred_cmip5(n,i)));
	end
	fprintf(fid,'\r\n');
end
fprintf(fid,'%04d-%04d',years_model(20),years_model(1));
for i=1:size(pred_cmip5,2)
	fprintf(fid,', %f',nanmean(exp(pred_cmip5(1:20,i))));
end
fprintf(fid,'\r\n');
fprintf(fid,'%04d-%04d',years_model(40),years_model(21));
for i=1:size(pred_cmip5,2)
	fprintf(fid,', %f',nanmean(exp(pred_cmip5(21:40,i))));
end
fprintf(fid,'\r\n');
fclose(fid);
%%%%%%%%%%%%
pred_lowice=pred_lowice+B1(4)*repmat(years_model(:),1,size(pred_lowice,2));
fid=fopen([workPath 'data/BA_lowice.csv'],'w');
fprintf(fid,'YYYY');
for i=1:size(pred_lowice,2)
	fprintf(fid,', GCM%02d-Low-Ice',i);
end
fprintf(fid,'\r\n');
for n=1:size(pred_lowice,1)
	fprintf(fid,'%04d',years_model(n));
	for i=1:size(pred_lowice,2)
		fprintf(fid,', %f',exp(pred_lowice(n,i)));
	end
	fprintf(fid,'\r\n');
end
fprintf(fid,'%04d-%04d',years_model(20),years_model(1));
for i=1:size(pred_lowice,2)
	fprintf(fid,', %f',nanmean(exp(pred_lowice(1:20,i))));
end
fprintf(fid,'\r\n');
fprintf(fid,'%04d-%04d',years_model(40),years_model(21));
for i=1:size(pred_lowice,2)
	fprintf(fid,', %f',nanmean(exp(pred_lowice(21:40,i))));
end
fprintf(fid,'\r\n');
fclose(fid);
%%%%%%%%%%%%
pred_ecearth=pred_ecearth+B1(4)*repmat(years_model(:),1,size(pred_ecearth,2));
fid=fopen([workPath 'data/BA_lowice_EC-EARTH.csv'],'w');
fprintf(fid,'YYYY');
for i=1:size(pred_ecearth,2)
	fprintf(fid,', GCM%02d-Low-Ice_EC-EARTH',i);
end
fprintf(fid,'\r\n');
for n=1:size(pred_ecearth,1)
	fprintf(fid,'%04d',years_model(n));
	for i=1:size(pred_ecearth,2)
		fprintf(fid,', %f',exp(pred_ecearth(n,i)));
	end
	fprintf(fid,'\r\n');
end
fprintf(fid,'%04d-%04d',years_model(20),years_model(1));
for i=1:size(pred_ecearth,2)
	fprintf(fid,', %f',nanmean(exp(pred_ecearth(1:20,i))));
end
fprintf(fid,'\r\n');
fprintf(fid,'%04d-%04d',years_model(40),years_model(21));
for i=1:size(pred_ecearth,2)
	fprintf(fid,', %f',nanmean(exp(pred_ecearth(21:40,i))));
end
fprintf(fid,'\r\n');
fclose(fid);
%%%%%%%%%%%%
fid=fopen([workPath 'data/BA_obs_model.csv'],'w');
fprintf(fid,'YYYY, Obs, Model \r\n');
for n=1:length(Y1)
	fprintf(fid,'%04d, %f, %f \r\n',years(n), exp(yor(n)), exp(Y1(n)));
end
fprintf(fid,'%04d-%04d, %f, %f \r\n',years(31),years(50), nanmean(exp(yor(31:50)')), nanmean(exp(Y1(31:50)')));
fprintf(fid,'%04d-%04d, %f, %f \r\n',years(51),years(70), nanmean(exp(yor(51:70)')), nanmean(exp(Y1(51:70)')));
fprintf(fid,'%04d-%04d, %f, %f \r\n',years(31),years(70), nanmean(exp(yor(31:70)')), nanmean(exp(Y1(31:70)')));
fclose(fid);
