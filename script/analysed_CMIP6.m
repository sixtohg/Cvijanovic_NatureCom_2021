clear all,clc,close all
%% cd('/home/sixto/Documents/MLToolbox_R2013/MeteoLab/'),init
workPath=[pwd '/../'];cd([workPath 'script/'])
dir_out='/Users/marco/Documents/output/ivana_ice_fires/';
domain='baileys';
dir_fire=[workPath 'data/'];
dir_sim_pr=[workPath 'data/CMIP5_CMIP6/pr/ssp585/'];
dir_sim_tasmax=[workPath 'data/CMIP5_CMIP6/tasmax/ssp585/'];

listGCMs_SPI = dir(fullfile(dir_sim_pr, 'pr_*_spi4.mat'));
listGCMs_TX = dir(fullfile(dir_sim_tasmax, 'tasmax*.mat'));

%% no std

SP_rcp85=zeros(40,length(listGCMs_SPI))*NaN;
TX_rcp85=zeros(40,length(listGCMs_TX))*NaN;


figure;

for ifile=1:length(listGCMs_SPI)
    %if isempty(strfind(listGCMs_SPI(ifile).name,'CMIP6_MIROC6'))
    
%     writeMat(file.path(dir_gcm, paste0('spi6_3_rcp45_',listGCMs[ifile],'_north_rcp45.mat')), spi6_3=spi6_3_rcp45)
%   writeMat(file.path(dir_gcm, paste0('tx_4_7_rcp45_',listGCMs[ifile],'_north_rcp45.mat')), tx_4_7=tx_4_7_rcp45)
%   writeMat(file.path(dir_gcm, paste0('spi6_3_rcp85_',listGCMs[ifile],'_north_rcp45.mat')), spi6_3=spi6_3_rcp85)
%   writeMat(file.path(dir_gcm, paste0('tx_4_7_rcp85_',listGCMs[ifile],'_north_rcp45.mat')), tx_4_7=tx_4_7_rcp85)
%   


%   

      load([dir_sim_pr listGCMs_SPI(ifile).name])
      SP_rcp85(:,ifile)=spi4_3;
      SP_changes_rcp85(ifile)=nanmean(spi4_3(21:end))-nanmean(spi4_3(1:20));
      
      load([dir_sim_tasmax listGCMs_TX(ifile).name])
      TX_rcp85(:,ifile)=tx_4_9;
      
      plot(tx_4_9)
      
      TX_changes_rcp85(ifile)=nanmean(tx_4_9(21:end))-nanmean(tx_4_9(1:20));
      
    %else
    %    figure;plot(tx_4_9)
    %    nanmean(tx_4_9(21:end))-nanmean(tx_4_9(1:20))
        
    %    figure;plot(spi4_3)
    %    nanmean(spi4_3(21:end))-nanmean(spi4_3(1:20))
    %end
end


figure; plot(TX_rcp85)
hold on;
listGCMs_TX(17).name
plot(TX_rcp85(:,17),'k-o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
cc



%% rcp85


savename=[dir_out 'SP_changes_rcp85_cmip6.mat']
save(savename,'SP_changes_rcp85')
savename=[dir_out 'TX_changes_rcp85_cmip6.mat']
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


nomeout=[dir_out,'boxplot_gcm_rcp85_cmip6'];
print( gcf, '-depsc2', nomeout ,'-painters')

figure;plot(SP_rcp85(:),TX_rcp85(:),'o')
%[a,b]=corr(SP_rcp85(:),TX_rcp85(:),'rows','complete')

