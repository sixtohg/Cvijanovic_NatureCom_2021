clear all
close all

%%

years=1951:2020;
dir_fire='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/';
dir_out='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/paper/figures_temp/';
domain='baileys';

%% frap 
savename = [dir_fire,'fire_1950_2020_',domain,'_mjjas_gt300.mat'];
%savename = [dir_fire,'fire_1950_2020_',domain,'_mjjaso.mat'];
load(savename,'FIRE')
anni_fire=1950:2020;
[C,IA,IB] = intersect(years,anni_fire);
FIRE=FIRE(IB);

%savename = [dir_fire,'fire_1950_2020_',domain,'_mjjaso_adjusted.mat'];
%load(savename,'FIRE')

%% mtbs

nomefile = [dir_fire,'fire_1984_2019_monthly_',domain,'.mat'];
load(nomefile,'BA')

anni_mtbs=1984:2019;
[C,IA,IB] = intersect(years,anni_mtbs);

BA=BA';
for iyear=1:length(anni_mtbs)
    i1=(iyear-1)*12+5
    i2=i1+4
    dum(iyear)=nansum(BA(i1:i2));
end

MTBS=zeros(length(years),1)*NaN;
MTBS(IA)=dum(IB); clear BAS

%% westerling

dir_west='/Users/marco/Documents/dati/fire_us/';
load([dir_west,'westerling_1980-2004_mjjas.mat']);
anni_west=1980:2004;
[C,IA,IB] = intersect(years,anni_west);
WEST=zeros(length(years),1)*NaN;
WEST(IA)=BAs(IB); clear BAs




%% 
FIRE=FIRE/1000;
MTBS=MTBS*10;
WEST=WEST/100000; 


%%
figure1 = figure
axes1 = axes('Parent',figure1);
plot(years,FIRE,'-o','color',[255,0,0]/255,'LineWidth',1.5)
hold on
plot(years,MTBS,'-s','color',[0,0,0]/255,'LineWidth',1.5)
plot(years,WEST*100,'-d','color',[100,100,100]/255,'LineWidth',1.5)
xlim([years(1)-1 years(end)+1])
legend('FRAP','MTBS','Westerling et al. (2003)','Location','NorthWest')
plot(years,FIRE,'-o','color',[255,0,0]/255,'LineWidth',1.5)
xlabel('Years','FontSize',20)
ylabel('Burned Area (1000 km^2)','FontSize',20)
set(gca,'FontSize',18)
xlim([1949 2022])
%set(gca, 'YScale', 'log')
%ylim([50 20000])


%% 
[r p]=corrcoef(FIRE,MTBS,'rows','complete')
[r p]=corrcoef(FIRE,WEST,'rows','complete')

file=[dir_out,'figure1a_temp.eps']
print( gcf, '-depsc2', file ,'-painters')
