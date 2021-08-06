clear all
close all

%%

years=1950:2020;
dir_fire='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/';
dir_out='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/paper/figures_temp/';
domain='baileys';

%% frap 1950-2020
nomefile = [dir_fire,'frap_baileys.mat'];
load(nomefile,'frap')
fires=frap;

min(fires(:,1))
max(fires(:,1))
ino=fires(:,1)>2020;
fires(ino,:)=[];

fires(:,4)=fires(:,4)*0.0040468564224; %acres to km2
[dum,iok]=sort(fires(:,1));

fires=fires(iok,:);
igt300=(fires(:,4)>=1.2141);
fires_gt300=fires(igt300,:);
%300 acres
%300*0.0040468564224=1.2141

figure;semilogy(fires(:,1),fires(:,4),'ko')
grid on
hold on
semilogy(xlim, [1.2141, 1.2141], 'Color', 'r', 'LineWidth', 2);
xlabel('Years','FontSize',18)
ylabel('Burned Area (km^2)','FontSize',18)
set(gca,'FontSize',18)
xlim([1949 2022])
file=[dir_out,'figure_gt300.eps']
print( gcf, '-depsc2', file ,'-painters')

%% monthly
% fire_2020[,1 ] = anni
% fire_2020[,2 ] = mesi
% fire_2020[,3 ] = giorni
% fire_2020[,4 ] = fires

FIREM=zeros(length(years),1);
FIREM_gt300=zeros(length(years),1);
k=0;
for iyear=1:length(years)
    for im=1:12
        k=k+1;
        iok=find(fires(:,1)==years(iyear) & fires(:,2)==im);
        FIREM(k)=sum(fires(iok,4));
        iok=find(fires_gt300(:,1)==years(iyear) & fires_gt300(:,2)==im);
        FIREM_gt300(k)=sum(fires_gt300(iok,4));
    end
end


figure;plot(FIREM)
hold on
plot(FIREM_gt300,'ro')



%% mjjas

for iyear=1:length(years)
    i1=(iyear-1)*12+5;
    i2=i1+4;
    FIRE2(iyear)=nansum(FIREM(i1:i2));
    FIRE3(iyear)=nansum(FIREM_gt300(i1:i2));
end




%savename = [dir_fire,'fire_1950_2020_',domain,'_mjjas.mat'];
%load(savename,'FIRE')


figure;plot(years,FIRE2/1000,'-o','color',[255,0,0]/255,'LineWidth',1.5)
hold on
plot(years,FIRE3/1000,'-s','color',[0,0,0]/255,'LineWidth',1.5)
%plot(FIRE,'ro')
xlabel('Years','FontSize',18)
ylabel('Burned Area (1000 km^2)','FontSize',18)
legend('All fires','Fires greater then 1.21 km^2','Location','NorthWest')
set(gca,'FontSize',18)
xlim([1949 2022])
file=[dir_out,'figure_gt300_mjjas.eps']
print( gcf, '-depsc2', file ,'-painters')

sum(FIRE3)/sum(FIRE2)

FIRE=FIRE3;
savename = [dir_fire,'fire_1950_2020_',domain,'_mjjas_gt300.mat'];
save(savename,'FIRE')