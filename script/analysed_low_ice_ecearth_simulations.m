clear all,clc,%% close all
%% cd('/home/sixto/Documents/MLToolbox_R2013/MeteoLab/'),init
workPath=[pwd '/../'];cd([workPath 'script/'])
dir_out='/Users/marco/Documents/output/ivana_ice_fires/';
domain='baileys';
dir_sim=[workPath 'data/low_ice_ecearth/'];

%% no std

SP=zeros(40,12)*NaN;
TX=zeros(40,12)*NaN;

for ifile=1:10
    
  %writeMat(file.path(dir_sim, paste0('spi3_4_',ifile-1,'_nostd_',domain,'.mat')), spi4_3=spi4_3)
  %writeMat(file.path(dir_sim, paste0('tx_4_10_',ifile-1,'_nostd_',domain,'.mat')), tx_4_10=tx_4_10)
  
      load([dir_sim 'spi4_3_' num2str(ifile-1) '_nostd_' char(domain) '.mat'])
      
      
      SP(:,ifile)=spi4_3;
      SP_changes_ecearth(ifile)=nanmean(spi4_3(21:end))-nanmean(spi4_3(1:20));
      
      load([dir_sim  'tx_4_9_' num2str(ifile-1) '_nostd_' char(domain) '.mat'])
      TX(:,ifile)=tx_4_9;
      TX_changes_ecearth(ifile)=nanmean(tx_4_9(21:end))-nanmean(tx_4_9(1:20))
end



%%

savename=[dir_out 'SP_changes_lowice_ecearth.mat']
save(savename,'SP_changes_ecearth')
savename=[dir_out 'TX_changes_lowice_ecearth.mat']
save(savename,'TX_changes_ecearth')


dataplot=zeros(1,size(SP_changes_ecearth,2),2)*NaN;


dataplot(1,1:size(SP_changes_ecearth,2),1)=SP_changes_ecearth'; 
dataplot(1,1:size(TX_changes_ecearth,2),2)=TX_changes_ecearth'; 


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


nomeout=[dir_out,'boxplot_ecearth_simulations'];
print( gcf, '-depsc2', nomeout ,'-painters')




figure;plot(SP(:),TX(:),'o')
[a,b]=corr(SP(:),TX(:),'rows','complete')


