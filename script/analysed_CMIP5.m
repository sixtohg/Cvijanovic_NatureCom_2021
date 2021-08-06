clear all,clc,%% close all
%% cd('/home/sixto/Documents/MLToolbox_R2013/MeteoLab/'),init
workPath=[pwd '/../'];cd([workPath 'script/'])
dir_out='/Users/marco/Documents/output/ivana_ice_fires/';
domain='baileys';
%% 
dir_sim_pr=[workPath 'data/CMIP5_CMIP6/pr/rcp85/'];
dir_sim_tasmax=[workPath 'data/CMIP5_CMIP6/tasmax/rcp85/'];
listGCMs={'ACCESS1-3','ACCESS1-3','BCC-CSM1-1','CNRM-CM5','CNRM-CM5','CanESM2','EC-EARTH','EC-EARTH','FGOALS-g2','GFDL-ESM2G', ...
    'GFDL-ESM2M','inmcm4','inmcm4','IPSL-CM5A-LR','MIROC-ESM','MIROC5','MPI-ESM-LR','MRI-CGCM3','NorESM1-M','NorESM1-M','HadGEM2-ES'};

SP_rcp85=zeros(40,length(listGCMs))*NaN;
TX_rcp85=zeros(40,length(listGCMs))*NaN;


for ifile=1:length(listGCMs)

      load([dir_sim_pr,'spi4_3_',num2str(ifile-1),'_nostd_',domain,'.mat'])
      SP_rcp85(:,ifile)=spi4_3;
      SP_changes_rcp85(ifile)=nanmean(spi4_3(21:end))-nanmean(spi4_3(1:20));
      
      load([dir_sim_tasmax,'tx_4_9_',num2str(ifile-1),'_nostd_',domain,'.mat'])
      TX_rcp85(:,ifile)=tx_4_9;
      TX_changes_rcp85(ifile)=nanmean(tx_4_9(21:end))-nanmean(tx_4_9(1:20))
end



%% rcp85


savename=[dir_out 'SP_changes_rcp85.mat']
save(savename,'SP_changes_rcp85')
savename=[dir_out 'TX_changes_rcp85.mat']
save(savename,'TX_changes_rcp85')


dataplot=zeros(1,size(SP_changes_rcp85,2),2)*NaN;

dataplot(1,1:size(SP_changes_rcp85,2),1)=SP_changes_rcp85'; 
dataplot(1,1:size(TX_changes_rcp85,2),2)=TX_changes_rcp85'; 


c=[150 150 150]/255;

hFig=figure;aboxplot3(dataplot,'colormap',c,'labels',{'SPI','TMAX'})
fntsz=20;

% for i=1:2
%     text(i,55,num2str(round(prctile(dataplot(1,:,i),50))),'FontSize',fntsz);
% end


ylabel('Ensemble mean differences','FontSize',fntsz);
xlabel('Stardardized climate variables','FontSize',fntsz);

gridxy([],[0:10:60],'Color','k','Linestyle',':');
set(gca,'FontSize',fntsz)


nomeout=[dir_out,'boxplot_gcm_rcp85_baileys_cal'];
print( gcf, '-depsc2', nomeout ,'-painters')

figure;plot(SP_rcp85(:),TX_rcp85(:),'o')
[a,b]=corr(SP_rcp85(:),TX_rcp85(:),'rows','complete')

