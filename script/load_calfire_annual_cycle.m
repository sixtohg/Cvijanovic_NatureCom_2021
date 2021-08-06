clear all
close all

%%

years=1950:2020;
dir_fire='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/';
dir_out='/Users/marco/Documents/output/ivana_ice_fires/';
domain='baileys';

%% frap
%fire18_1_JJA_1950_2019_baileys_cal.mat
%fire_1950_2019_monthly_baileys.mat
nomefile = [dir_fire,'fire_1950_2020_monthly_',domain,'.mat'];
load(nomefile,'BA')
BA_1950_2020=BA*0.0040468564224/1000;
clear BA

nomefile = [dir_fire,'fire_1950_2020_monthly_cal.mat'];
load(nomefile,'BA')
BA_1950_2020_cal=BA*0.0040468564224/1000;
clear BA
BAm_cal=BA_1950_2020_cal;


%%
BAm=BA_1950_2020;



savename = [dir_fire,'fire_1950_2020_',domain,'_monthly.mat'];
save(savename,'BAm')




for im=1:12
    BAmm(im,:)=(BAm(im:12:end));
end




%% --> ba concentrata nei mesi tra 6 e 10
%ba in jjaso
nansum(nansum(BAmm(5:10,:),1),2)/nansum(nansum(BAmm(:,:),1),2)
%ba in jja
nansum(nansum(BAmm(6:8,:),1),2)/nansum(nansum(BAmm(:,:),1),2)
nansum(nansum(BAmm(6:9,:),1),2)/nansum(nansum(BAmm(:,:),1),2)
nansum(nansum(BAmm(5:9,:),1),2)/nansum(nansum(BAmm(:,:),1),2)


color_fire=[255,69,0]/255;


figure;
aboxplot2(BAmm','colormap',color_fire,'labels',[1:12]); % Advanced box plot
xlabel('Months','FontSize',20)
ylabel('Burned Area (1000 km^2)','FontSize',20)
set(gca,'FontSize',18)
title('BA','FontSize',24)
%ylim([0 4100])
file=[dir_out,'annual_cycle_FRAP.eps']
print( gcf, '-depsc2', file ,'-painters')


%% save fire series 1950-2020

BAm=BAm';
for iyear=1:length(years)
    i1=(iyear-1)*12+5
    i2=i1+5
    FIRE(iyear)=nansum(BAm(i1:i2));
end


figure;plot(FIRE)

savename = [dir_fire,'fire_1950_2020_',domain,'_mjjaso.mat'];
save(savename,'FIRE')



%% annual bal vs cal
for iyear=1:length(years)
    i1=(iyear-1)*12+1;
    i2=i1+11;
    FIRE(iyear)=nansum(BAm(i1:i2));
end

BAm_cal=BAm_cal';
for iyear=1:length(years)
    i1=(iyear-1)*12+1;
    i2=i1+11;
    FIRE_cal(iyear)=nansum(BAm_cal(i1:i2));
end


%FIRE=FIRE(end-39:end);
%FIRE_cal=FIRE_cal(end-39:end);
figure; plot(FIRE); hold on; plot(FIRE_cal,'r')

figure; plot(FIRE./FIRE_cal)

sum(FIRE)/sum(FIRE_cal)


%% mjjas

for iyear=1:length(years)
    i1=(iyear-1)*12+5;
    i2=i1+4;
    FIRE(iyear)=nansum(BAm(i1:i2));
end


figure;plot(FIRE)

savename = [dir_fire,'fire_1950_2020_',domain,'_mjjas.mat'];
save(savename,'FIRE')

figure;plot(FIRE(end-39:end))
figure;plot(log(FIRE(end-39:end)))