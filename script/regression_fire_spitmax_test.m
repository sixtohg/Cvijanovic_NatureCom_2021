clear all
close all

%%

years=(1951:2020);
alpha = 0.05;
NB=1000;
dir_fire='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/';
dir_out='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/paper/figures_temp/';
domain='baileys';

%%
savename = [dir_fire,'fire_1950_2020_',domain,'_mjjas_gt300.mat'];
load(savename,'FIRE')
anni_fire=1950:2020;
[C,IA,IB] = intersect(years,anni_fire);
FIRE=FIRE(IB)/1000;




figure; plot(FIRE)


yor=(log(FIRE)');

figure; plot(yor)




%% califonia nclimgrid
filename = [dir_fire,'spi4_3_nclimgrid.mat'];
load(filename);
filename = [dir_fire,'tx_4_9_nclimgrid.mat'];
load(filename);


SP=scale_1981_2000(spi4_3,years);
TX=scale_1981_2000(tx_4_9,years);

[a,b]=corrcoef(SP,TX)



%predit={ '$SPI4_{3}$ & '  '$Tx_{(4-10)}$ & ' };
%nomepredit={ 'spi4_3' 'tx_4_10'  }';

%% only spi

Xorig1 = [ones(size(SP,1),1) nandetrend(SP)  years'  ]; %osservati
[B1,BINT,R,RINT,STATS1] = regress((yor),Xorig1,alpha);
STATS1
[~,akc]=aic2(length(find(isnan(yor)==0)),size(Xorig1,2),nanvar(R),nanvar((yor)))


%% only TSMAX

Xorig1 = [ones(size(SP,1),1) nandetrend(TX)  years'  ]; %osservati
[B1,BINT,R,RINT,STATS1] = regress((yor),Xorig1,alpha);
STATS1
[~,akc]=aic2(length(find(isnan(yor)==0)),size(Xorig1,2),nanvar(R),nanvar((yor)))

%% regression detrended climate with t^2

Xorig1 = [ones(size(SP,1),1) (SP) (TX) years'  ]; %osservati
[B1,BINT,R,RINT,STATS1] = regress((yor),Xorig1,alpha);
STATS1
[~,akc]=aic2(length(find(isnan(yor)==0)),size(Xorig1,2),nanvar(R),nanvar((yor)))

Xorig1 = [ones(size(SP,1),1) (SP) (TX) years' (years').^2 (years').^3]; %osservati
[B1,BINT,R,RINT,STATS1] = regress((yor),Xorig1,alpha);
STATS1
[~,akc]=aic2(length(find(isnan(yor)==0)),size(Xorig1,2),nanvar(R),nanvar((yor)))


%Xorig1 = [ones(size(SP,1),1) nandetrend(SP) nandetrend(TX) scale_1981_2000(years',years)  ]; %osservati
%Xorig1 = [ones(size(SP,1),1) nandetrend(SP) nandetrend(TX,2) years' (years').^2 ]; %osservati
Xorig1 = [ones(size(SP,1),1) (SP) (TX) years' (years').^2 ]; %osservati
[B1,BINT,R,RINT,STATS1] = regress((yor),Xorig1,alpha);
STATS1
[~,akc]=aic2(length(find(isnan(yor)==0)),size(Xorig1,2),nanvar(R),nanvar((yor)))

Y1 = (Xorig1*B1);
%Y2 = B1(1) + B1(2) * SP + B1(3) * TX + B1(4) * scale_1981_2000(years',years);

[a,b]=corrcoef((yor),Y1)
[a,b]=corrcoef(nandetrend(SP),nandetrend(TX))


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
format long g
B1

100*(exp(B1(3))-1)
100*(exp(-B1(2))-1)
cc

%% regression detrended climate



%Xorig1 = [ones(size(SP,1),1) nandetrend(SP) nandetrend(TX) scale_1981_2000(years',years)  ]; %osservati
Xorig1 = [ones(size(SP,1),1) nandetrend(SP) nandetrend(TX) years'  ]; %osservati
[B1,BINT,R,RINT,STATS1] = regress((yor),Xorig1,alpha);
STATS1
[~,akc]=aic2(length(find(isnan(yor)==0)),size(Xorig1,2),nanvar(R),nanvar((yor)))

Y1 = (Xorig1*B1);
%Y2 = B1(1) + B1(2) * SP + B1(3) * TX + B1(4) * scale_1981_2000(years',years);

[a,b]=corrcoef((yor),Y1)
[a,b]=corrcoef(nandetrend(SP),nandetrend(TX))


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


B1

B1(3)

100*(exp(B1(3))-1)
cc

%% changes



%% regression all drivers and only climate drivers


Xorig1 = [ones(size(SP,1),1) SP TX scale_1981_2000(years',years)  ]; %osservati
[B1,BINT,R,RINT,STATS1] = regress((yor),Xorig1,alpha);
STATS1
[~,akc]=aic2(length(find(isnan(yor)==0)),size(Xorig1,2),nanvar(R),nanvar(nanpstd(yor)));

Y1 = (Xorig1*B1);
Y2 = B1(1) + B1(2) * SP + B1(3) * TX + B1(4) * scale_1981_2000(years',years);

[a,b]=corrcoef((yor),Y1)


figure;
plot(years,exp(yor),'-ko','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
hold on
plot(years,exp(Y1),'k:o','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',6)
plot(years,exp(Y2),'r:o','LineWidth',1,'MarkerEdgeColor','r','MarkerSize',6)
legend('Observed','Modelled','Location','Best')
xlabel('Years') 

(nanmean(Y1)-nanmean(yor))/nanmean(yor)
(nanmean(Y1)-nanmean(yor))
nanmean(yor)
nanmean(Y1)
nanmean(exp(yor))
nanmean(FIRE,2)
exp(nanmean(yor))
exp(nanmean(Y1))
100*(nanmean(exp(Y1))-nanmean(exp(yor)))/nanmean(exp(yor))




bootdat = [(yor),Xorig1];
bootb1 = bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat);
bootb1 = sort(bootb1);
bootCI1 = [bootb1(25,:); bootb1(975,:)]  %Report 5% confidence intervals
B1

clear bootdat

for ib=1:1000
    all_drivers(:,ib) = bootb1(ib,1) + bootb1(ib,2) * SP + bootb1(ib,3)  * TX + bootb1(ib,4) * scale_1981_2000(years',years);
    climate_drivers(:,ib) = bootb1(ib,1) + bootb1(ib,2) * SP + bootb1(ib,3)  * TX;    
end


data(:,1)=nanmedian(all_drivers,2);
data(:,2)=prctile(all_drivers',2.5); %min(YRall');
data(:,3)=prctile(all_drivers',97.5); %max(YRall');

%dataR(:,1)=mean(Ycall,2);
dataR(:,1)=nanmedian(climate_drivers,2);
dataR(:,2)=prctile(climate_drivers',2.5); %min(YRall');
dataR(:,3)=prctile(climate_drivers',97.5); %max(YRall');


colorsg    = [1,0,0]; 
colorsgE    = [255, 165, 0]/255; %grey
widths= [2 0 0]';


figure;
drawSpread(exp(data),'xvalues',years','lines','no','indexes',(1:3)','widths',widths,'colorsg',colorsgE);
hold on
drawSpread(exp(dataR),'xvalues',years','lines','no','indexes',(1:3)','widths',widths,'colorsg',colorsg);
plot(years,exp(yor),'-ko','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)

legend('All drivers','Climate drivers','Observed','Location','Best')
xlabel('Years','FontSize',18) 
ylabel('log(BA) [log(1000 km^2)]','FontSize',18) 
%gridxy(years(ith),'Color','k','Linestyle',':') ;
set(gca,'FontSize',18,'xlim',[years(1)-1 years(end)+1])
set(gcf, 'PaperPositionMode', 'auto','renderer', 'painters');
file=[dir_out,'figure1x_temp.eps']
print( gcf, '-depsc2', file ,'-painters')






%% barplot

% period1=1951:2020;
% [C,IA,IB] = intersect(period1,years);
% FIRE1=mean(FIRE(IB))
% ALLD1=mean(exp(data(IB,1)))
% ALLD1up=mean(exp(data(IB,3)))
% ALLD1low=mean(exp(data(IB,2)))
% CLID1=mean(exp(dataR(IB,1)))
% CLID1up=mean(exp(dataR(IB,3)))
% CLID1low=mean(exp(dataR(IB,2)))

% 
% period2=1981:2000;
% [C,IA,IB] = intersect(period2,years);
% FIRE2=mean(FIRE(IB));
% ALLD2=mean(exp(data(IB,1)))
% ALLD2up=mean(exp(data(IB,3)))
% ALLD2low=mean(exp(data(IB,2)))
% CLID2=mean(exp(dataR(IB,1)))
% CLID2up=mean(exp(dataR(IB,3)))
% CLID2low=mean(exp(dataR(IB,2)))
% 
% period3=2001:2020;
% [C,IA,IB] = intersect(period3,years);
% FIRE3=mean(FIRE(IB));
% ALLD3=mean(exp(data(IB,1)))
% ALLD3up=mean(exp(data(IB,3)))
% ALLD3low=mean(exp(data(IB,2)))
% CLID3=mean(exp(dataR(IB,1)))
% CLID3up=mean(exp(dataR(IB,3)))
% CLID3low=mean(exp(dataR(IB,2)))

%%

%y=[FIRE1 ALLD1 CLID1; FIRE2 ALLD2 CLID2; FIRE3 ALLD3 CLID3];
%y=[FIRE1 ALLD1 CLID1];
%err=ones(1,3,2)*NaN 
%err(1,:,1) = [0 ALLD1-ALLD1low CLID1-CLID1low];
%err(1,:,2) = [0 ALLD1up-ALLD1 CLID1up-CLID1];
% err(2,:,1) = [0 ALLD2-ALLD2low CLID2-CLID2low];
% err(2,:,2) = [0 ALLD2up-ALLD2 CLID2up-CLID2];
% err(3,:,1) = [0 ALLD3-ALLD3low CLID3-CLID3low];
% err(3,:,2) = [0 ALLD3up-ALLD3 CLID3up-CLID3];


% figure(1); clf; 
% hb = bar(y,'FaceColor','flat'); % get the bar handles
% 
% legend('Observed','All drivers','Climate drivers','Location','Best')
% hold on;
% for k = 1:size(y,2)
%     
%     %hb(k).CData = k;
%     
%     % get x positions per group
%     xpos = hb(k).XData + hb(k).XOffset;
%     
%     
%     % draw errorbar
%     errorbar(xpos, y(:,k), err(:,k,1),err(:,k,2), 'LineStyle', 'none', ... 
%         'Color', 'k', 'LineWidth', 1);
% end

% % % Set Axis properties
% %set(gca,'xticklabel',{'1951-2000'; '1981-2000'; '20001-2020'});
% set(gca,'xticklabel',{'1951-2020'});
% % ylim([200 360])
% ylabel('mean BA/ year [1000 km^2]')
% xlabel('periods')
% set(gca,'FontSize',18)
% set(gcf, 'PaperPositionMode', 'auto','renderer', 'painters');
% file=[dir_out,'figure1y_temp.eps']
% print( gcf, '-depsc2', file ,'-painters')

%% barplot
data(:,1)=nanmedian(all_drivers,2);
data(:,2)=prctile(all_drivers',2.5); %min(YRall');
data(:,3)=prctile(all_drivers',97.5); %max(YRall');

dataR(:,1)=nanmedian(climate_drivers,2);
dataR(:,2)=prctile(climate_drivers',2.5); %min(YRall');
dataR(:,3)=prctile(climate_drivers',97.5); %max(YRall');

period1=1951:2000;
[C,IA,IB] = intersect(period1,years);
FIRE1=mean(FIRE(IB));
ALLD1=mean(exp(data(IB,1)));
ALLD1up=mean(exp(data(IB,3)));
ALLD1low=mean(exp(data(IB,2)));
CLID1=mean(exp(dataR(IB,1)));
CLID1up=mean(exp(dataR(IB,3)));
CLID1low=mean(exp(dataR(IB,2)));

period2=1981:2000;
[C,IA,IB] = intersect(period2,years);
FIRE2=mean(FIRE(IB));
ALLD2=mean(exp(data(IB,1)));
ALLD2up=mean(exp(data(IB,3)));
ALLD2low=mean(exp(data(IB,2)));
CLID2=mean(exp(dataR(IB,1)));
CLID2up=mean(exp(dataR(IB,3)));
CLID2low=mean(exp(dataR(IB,2)));

period3=2001:2020;
[C,IA,IB] = intersect(period3,years);
FIRE3=mean(FIRE(IB));
ALLD3=mean(exp(data(IB,1)));
ALLD3up=mean(exp(data(IB,3)));
ALLD3low=mean(exp(data(IB,2)));
CLID3=mean(exp(dataR(IB,1)));
CLID3up=mean(exp(dataR(IB,3)));
CLID3low=mean(exp(dataR(IB,2)));

y=[FIRE1 ALLD1 CLID1; FIRE2 ALLD2 CLID2; FIRE3 ALLD3 CLID3];
err=ones(2,3,2)*NaN ;
err(1,:,1) = [0 ALLD1-ALLD1low CLID1-CLID1low];
err(1,:,2) = [0 ALLD1up-ALLD1 CLID1up-CLID1];
err(2,:,1) = [0 ALLD2-ALLD2low CLID2-CLID2low];
err(2,:,2) = [0 ALLD2up-ALLD2 CLID2up-CLID2];
err(3,:,1) = [0 ALLD3-ALLD3low CLID3-CLID3low];
err(3,:,2) = [0 ALLD3up-ALLD3 CLID3up-CLID3];

figure(1); clf; 
hb = bar(y,'FaceColor','flat'); % get the bar handles

hb(1).CData(:,2)=0;
hb(1).CData(:,3)=1;
hb(1).CData(:,3)=0;

legend('Observed','All drivers','Climate drivers','Location','Best')
hold on;
for k = 1:size(y,2)
    %% hb(k).CData = k;
    %% get x positions per group
    xpos = hb(k).XData + hb(k).XOffset;
    %% draw errorbar
    errorbar(xpos, y(:,k), err(:,k,1),err(:,k,2), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);
end

%% Set Axis properties
set(gca,'xticklabel',{'1951-2000'; '1981-2000'; '2001-2020'});
% ylim([200 360])
ylabel('mean BA/ year [1000 km^2]')
xlabel('periods')
set(gca,'FontSize',18)
set(gcf, 'PaperPositionMode', 'auto','renderer', 'painters');
file=[dir_out,'figure1y_temp.eps']
print( gcf, '-depsc2', file ,'-painters')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% echale un vistazo a este c??digo:
%% 
%% regression_fire_spitmax_test.m
%% 
%% revisa que tambi??n el segunda gr??fico es correcto por favor, y si puede cambia los colores de las barras:
%% 
%% obs dir??a negro
%% y la otras dos, como la grafica de antes,  naranja y rojo, o como veas tu mejor
%% y a parte de todo, se entiende que quiere decir esta analisis? tu crees que vale la pena a??adirla?
%% 
%% gracias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all, clc, close all, fclose all;
cd('/home/sixto/Documents/MLToolbox_R2013/MeteoLab/'),init
workPath='/home/sixto/Dropbox/ivana_ice_fires/review/';cd([workPath 'script/'])
years=(1951:2020);
alpha = 0.05;
NB=1000;
dir_fire=[workPath 'data/'];
dir_out=[workPath 'paper/figures_temp/'];
%% dir_fire='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/';
%% dir_out='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/paper/figures_temp/';
domain='baileys';

savename = [dir_fire,'fire_1950_2020_',domain,'_mjjas_gt300.mat'];
load(savename,'FIRE')
anni_fire=1950:2020;
[C,IA,IB] = intersect(years,anni_fire);
FIRE=FIRE(IB)/1000;

figure; plot(years(IA),FIRE)

yor=(log(FIRE)');
figure; plot(years(IA),yor)

%% califonia nclimgrid
filename = [dir_fire,'spi4_3_nclimgrid.mat'];
load(filename);
filename = [dir_fire,'tx_4_9_nclimgrid.mat'];
load(filename);

SP=scale_1981_2000(spi4_3,years);
TX=scale_1981_2000(tx_4_9,years);

[a,b]=corrcoef(SP,TX);

%% predit={ '$SPI4_{3}$ & '  '$Tx_{(4-10)}$ & ' };
%% nomepredit={ 'spi4_3' 'tx_4_10'  }';

%% regression all drivers and only climate drivers
Xorig1 = [ones(size(SP,1),1) SP TX scale_1981_2000(years',years)]; %osservati
[B1,BINT,R,RINT,STATS1] = regress((yor),Xorig1,alpha);
%% STATS1
%% [~,akc]=aic2(length(find(isnan(yor)==0)),size(Xorig1,2),nanvar(R),nanvar(nanpstd(yor)));

Y1 = (Xorig1*B1);%% Y2 is the same but extended: Y2 = B1(1) + B1(2) * SP + B1(3) * TX + B1(4) * scale_1981_2000(years',years);
Y2 = B1(1) + B1(2) * SP + B1(3) * TX;%% Y2 including only the climate drivers
[a1,b1]=corrcoef((yor),Y1);
[a2,b2]=corrcoef((yor),Y2);

figure;
plot(years,(yor),'-ko','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
hold on
plot(years,(Y1),'r:o','LineWidth',1,'MarkerEdgeColor','r','MarkerSize',6)
plot(years,(Y2),'b:o','LineWidth',1,'MarkerEdgeColor','b','MarkerSize',6)
legend('Observed','Modelled (All drivers)','Modelled (Climate drivers)','Location','Best')
xlabel('Years') 

bootdat = [(yor),Xorig1];
bootb1 = bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat);
bootb1 = sort(bootb1);
bootCI1 = [bootb1(25,:); bootb1(975,:)];% Report 5% confidence intervals
%% B1

clear bootdat

for ib=1:1000
    all_drivers(:,ib) = bootb1(ib,1) + bootb1(ib,2) * SP + bootb1(ib,3)  * TX + bootb1(ib,4) * scale_1981_2000(years',years);
    climate_drivers(:,ib) = bootb1(ib,1) + bootb1(ib,2) * SP + bootb1(ib,3)  * TX;    
end
data=[all_drivers climate_drivers];
indices{1}=[1:1000];%% All drivers
indices{2}=[1001:2000];%% climate drivers
boundary=cell(2,1);boundary{1}={'prc2.5','prc97.5'}; boundary{2}={'prc2.5','prc97.5'};
%% h=figure;drawSpread(data,'xvalues',years','indexesg',indices,'colorsg',[255 0 0;101 101 101]/255,'lines','no','boundary',boundary,'alphasg',[0.5 0.5]); hold on
h=figure;drawSpread(data,'xvalues',years','indexesg',indices,'colorsg',[255 0 0;255, 165, 0]/255,'lines','no','widths',[2 0 0;2 0 0],'boundary',boundary,'alphasg',[0.5 0.5]); hold on
plot(years,(yor),'-ko','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)

legend('All drivers','Climate drivers','Observed','Location','Best')
xlabel('Years','FontSize',18) 
ylabel('log(BA) [log(1000 km^2)]','FontSize',18) 
%gridxy(years(ith),'Color','k','Linestyle',':') ;
set(gca,'FontSize',18,'xlim',[years(1)-1 years(end)+1])
set(gcf, 'PaperPositionMode', 'auto','renderer', 'painters');
file=[dir_out,'figure1x_temp.eps'];
print( gcf, '-depsc2', file ,'-painters')


%% barplot
data(:,1)=nanmedian(all_drivers,2);
data(:,2)=prctile(all_drivers',2.5); %min(YRall');
data(:,3)=prctile(all_drivers',97.5); %max(YRall');

dataR(:,1)=nanmedian(climate_drivers,2);
dataR(:,2)=prctile(climate_drivers',2.5); %min(YRall');
dataR(:,3)=prctile(climate_drivers',97.5); %max(YRall');

period1=1951:2000;
[C,IA,IB] = intersect(period1,years);
FIRE1=mean(FIRE(IB));
ALLD1=mean(exp(data(IB,1)));
ALLD1up=mean(exp(data(IB,3)));
ALLD1low=mean(exp(data(IB,2)));
CLID1=mean(exp(dataR(IB,1)));
CLID1up=mean(exp(dataR(IB,3)));
CLID1low=mean(exp(dataR(IB,2)));

period2=1981:2000;
[C,IA,IB] = intersect(period2,years);
FIRE2=mean(FIRE(IB));
ALLD2=mean(exp(data(IB,1)));
ALLD2up=mean(exp(data(IB,3)));
ALLD2low=mean(exp(data(IB,2)));
CLID2=mean(exp(dataR(IB,1)));
CLID2up=mean(exp(dataR(IB,3)));
CLID2low=mean(exp(dataR(IB,2)));

period3=2001:2020;
[C,IA,IB] = intersect(period3,years);
FIRE3=mean(FIRE(IB));
ALLD3=mean(exp(data(IB,1)));
ALLD3up=mean(exp(data(IB,3)));
ALLD3low=mean(exp(data(IB,2)));
CLID3=mean(exp(dataR(IB,1)));
CLID3up=mean(exp(dataR(IB,3)));
CLID3low=mean(exp(dataR(IB,2)));

y=[FIRE1 ALLD1 CLID1; FIRE2 ALLD2 CLID2; FIRE3 ALLD3 CLID3];
err=ones(2,3,2)*NaN ;
err(1,:,1) = [0 ALLD1-ALLD1low CLID1-CLID1low];
err(1,:,2) = [0 ALLD1up-ALLD1 CLID1up-CLID1];
err(2,:,1) = [0 ALLD2-ALLD2low CLID2-CLID2low];
err(2,:,2) = [0 ALLD2up-ALLD2 CLID2up-CLID2];
err(3,:,1) = [0 ALLD3-ALLD3low CLID3-CLID3low];
err(3,:,2) = [0 ALLD3up-ALLD3 CLID3up-CLID3];

figure(1); clf; 
hb = bar(y,'FaceColor','flat'); % get the bar handles

hb(1).CData(:,2)=0;
hb(1).CData(:,3)=1;
hb(1).CData(:,3)=0;

legend('Observed','All drivers','Climate drivers','Location','Best')
hold on;
for k = 1:size(y,2)
    %% hb(k).CData = k;
    %% get x positions per group
    xpos = hb(k).XData + hb(k).XOffset;
    %% draw errorbar
    errorbar(xpos, y(:,k), err(:,k,1),err(:,k,2), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);
end

%% Set Axis properties
set(gca,'xticklabel',{'1951-2000'; '1981-2000'; '2001-2020'});
% ylim([200 360])
ylabel('mean BA/ year [1000 km^2]')
xlabel('periods')
set(gca,'FontSize',18)
set(gcf, 'PaperPositionMode', 'auto','renderer', 'painters');
file=[dir_out,'figure1y_temp.eps']
print( gcf, '-depsc2', file ,'-painters')

%% Te pido otra cosa, que seguro que lo puedes hacer rapido,
clear all, clc, close all, fclose all;
cd('/home/sixto/Documents/MLToolbox_R2013/MeteoLab/'),init
workPath='/home/sixto/Dropbox/ivana_ice_fires/review/';cd([workPath 'script/'])
years=(1951:2020);
alpha = 0.05;
NB=1000;
dir_fire=[workPath 'data/'];
dir_out=[workPath 'paper/figures_temp/'];
%% dir_fire='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/';
%% dir_out='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/paper/figures_temp/';
domain='baileys';

savename = [dir_fire,'fire_1950_2020_',domain,'_mjjas_gt300.mat'];
load(savename,'FIRE')
anni_fire=1950:2020;
[C,IA,IB] = intersect(years,anni_fire);
FIRE=FIRE(IB)/1000;

%% figure; plot(years(IA),FIRE)
cFIRE= cumsum(FIRE');
%% valor cumumulado hasta el 2000 desde el 1951, obs
%% valor cumumulado hasta el 2000 desde el 1981, obs
%% valor cumumulado hasta el 2020 desde el 2001, obs
%% tabla=[nansum(FIRE(1:50)') nansum(FIRE(31:50)') nansum(FIRE(51:70)')];
tabla=[cFIRE(50) cFIRE(50)-cFIRE(30) cFIRE(70)-cFIRE(50)];

%% valor cumumulado hasta el 2050 desde el 2031, con los 4 grupos de modelos y intervalo de confianza al 95%, sin bootstrap
clear all, clc, close all, fclose all;
cd('/home/sixto/Documents/MLToolbox_R2013/MeteoLab/'),init
workPath='/home/sixto/Dropbox/ivana_ice_fires/review/';cd([workPath 'script/'])
years=(1951:2020);
alpha = 0.05;
NB=1000;
dir_fire=[workPath 'data/'];
dir_out=[workPath 'paper/figures_temp/'];
%% dir_fire='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/';
%% dir_out='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/paper/figures_temp/';
domain='baileys';

savename = [dir_fire,'fire_1950_2020_',domain,'_mjjas_gt300.mat'];
load(savename,'FIRE')
anni_fire=1950:2020;
[C,IA,IB] = intersect(years,anni_fire);
FIRE=FIRE(IB)/1000;
yor=(log(FIRE)');

%% califonia nclimgrid
filename = [dir_fire,'spi4_3_nclimgrid.mat'];load(filename);
filename = [dir_fire,'tx_4_9_nclimgrid.mat'];load(filename);
SP=scale_1981_2000(spi4_3,years);
TX=scale_1981_2000(tx_4_9,years);

%% regression all drivers and only climate drivers
Xorig1 = [ones(size(SP,1),1) SP TX scale_1981_2000(years',years)]; %osservati
[B1,BINT,R,RINT,STATS1] = regress((yor),Xorig1,alpha);

%% CMIP5:
dir_sim_pr='/home/sixto/Documents/Publications/Articulos/2020_Turco/scripts/GB-models/CMIP5_CMIP6/CMIP5_CMIP6/pr/rcp85/';
dir_sim_tasmax='/home/sixto/Documents/Publications/Articulos/2020_Turco/data/tasmax/rcp85/';
listGCMs_SPI = dir(fullfile(dir_sim_pr, 'pr_*.mat'));
listGCMs_TX = dir(fullfile(dir_sim_tasmax, 'tasmax*_max.mat'));
indGCMs=[2 2 3 11 11 22 28 28 33 34 35 59 59 41 48 49 52 56 57 57 37];
listGCMs_SPI=listGCMs_SPI(indGCMs);
listGCMs_TX=listGCMs_TX([1 1 15 2 2 3 4 4 5 6 7 16 16 9 10 11 12 13 14 14 8]);
refTmax=repmat(NaN,20,length(listGCMs_SPI));
refSPI=repmat(NaN,20,length(listGCMs_SPI));
for ifile=1:length(listGCMs_SPI)
	load(sprintf('/home/sixto/Dropbox/ivana_ice_fires/review/data/CMIP5/spi4_3_%d_nostd_baileys.mat',ifile-1));
	load(sprintf('/home/sixto/Dropbox/ivana_ice_fires/review/data/CMIP5/tx_4_9_%d_nostd_baileys.mat',ifile-1));
    tx_4_9(find(isnan(tx_4_9.*spi4_3)))=NaN;spi4_3(find(isnan(tx_4_9.*spi4_3)))=NaN;
	spi4_3=scale_1981_2000(spi4_3,[1981:2000 2031:2050]);
	tx_4_9=scale_1981_2000(tx_4_9,[1981:2000 2031:2050]);
	refTmax(:,ifile)=tx_4_9(1:20);refSPI(:,ifile)=spi4_3(1:20);
end
ba_cmip5_hist=exp(B1(1) + B1(2)*refSPI + B1(3)*refTmax);
ba_cmip5_hist_all=exp(B1(1) + B1(2)*refSPI + B1(3)*refTmax + B1(4) * scale_1981_2000([1981:2000]',[1981:2000]));

%% CMIP5:
dir_sim_pr='/home/sixto/Documents/Publications/Articulos/2020_Turco/scripts/GB-models/CMIP5_CMIP6/CMIP5_CMIP6/pr/rcp85/';
dir_sim_tasmax='/home/sixto/Documents/Publications/Articulos/2020_Turco/data/tasmax/rcp85/';
listGCMs_SPI = dir(fullfile(dir_sim_pr, 'pr_*.mat'));
listGCMs_TX = dir(fullfile(dir_sim_tasmax, 'tasmax*_max.mat'));
indGCMs=[2 2 3 11 11 22 28 28 33 34 35 59 59 41 48 49 52 56 57 57 37];
listGCMs_SPI=listGCMs_SPI(indGCMs);
listGCMs_TX=listGCMs_TX([1 1 15 2 2 3 4 4 5 6 7 16 16 9 10 11 12 13 14 14 8]);
ba_cmip5_proj=repmat(NaN,20,length(listGCMs_SPI));
for ifile=1:length(listGCMs_SPI)
	load(sprintf('/home/sixto/Dropbox/ivana_ice_fires/review/data/CMIP5/spi4_3_%d_nostd_baileys.mat',ifile-1));
	load(sprintf('/home/sixto/Dropbox/ivana_ice_fires/review/data/CMIP5/tx_4_9_%d_nostd_baileys.mat',ifile-1));
    tx_4_9(find(isnan(tx_4_9.*spi4_3)))=NaN;spi4_3(find(isnan(tx_4_9.*spi4_3)))=NaN;
	spi4_3=scale_1981_2000(spi4_3,[1981:2000 2031:2050]);
	tx_4_9=scale_1981_2000(tx_4_9,[1981:2000 2031:2050]);
	ba_cmip5_proj(:,ifile)=exp(B1(1) + B1(2)*spi4_3(21:40) + B1(3)*tx_4_9(21:40));
end

%% CMIP6:
dir_sim_pr='/home/sixto/Documents/Publications/Articulos/2020_Turco/scripts/GB-models/CMIP5_CMIP6/CMIP5_CMIP6/pr/ssp585/';
dir_sim_tasmax='/home/sixto/Documents/Publications/Articulos/2020_Turco/data/tasmax/ssp585/';
listGCMs_SPI = dir(fullfile(dir_sim_pr, 'pr_*.mat'));
listGCMs_TX = dir(fullfile(dir_sim_tasmax, 'tasmax*_max.mat'));
ba_cmip6_proj=repmat(NaN,20,length(listGCMs_SPI));
for ifile=1:length(listGCMs_SPI)
	load(sprintf('/home/sixto/Dropbox/ivana_ice_fires/review/data/CMIP6/spi4_3_%d_nostd_baileys.mat',ifile-1));
	load(sprintf('/home/sixto/Dropbox/ivana_ice_fires/review/data/CMIP6/tx_4_9_%d_nostd_baileys.mat',ifile-1));
    tx_4_9(find(isnan(tx_4_9.*spi4_3)))=NaN;spi4_3(find(isnan(tx_4_9.*spi4_3)))=NaN;
	spi4_3=scale_1981_2000(spi4_3,[1981:2000 2031:2050]);
	tx_4_9=scale_1981_2000(tx_4_9,[1981:2000 2031:2050]);
	ba_cmip6_proj(:,ifile)=exp(B1(1) + B1(2)*spi4_3(21:40) + B1(3)*tx_4_9(21:40));
end

%% Low-Ice:
ba_low_ice_proj=repmat(NaN,20,12);
for ifile=1:12
	load(sprintf('/home/sixto/Dropbox/ivana_ice_fires/review/data/low_ice/spi4_3_%d_nostd_baileys.mat',ifile-1));
	load(sprintf('/home/sixto/Dropbox/ivana_ice_fires/review/data/low_ice/tx_4_9_%d_nostd_baileys.mat',ifile-1));
    tx_4_9(find(isnan(tx_4_9.*spi4_3)))=NaN;spi4_3(find(isnan(tx_4_9.*spi4_3)))=NaN;
	spi4_3=scale_1981_2000(spi4_3,[1981:2000 2031:2050]);
	tx_4_9=scale_1981_2000(tx_4_9,[1981:2000 2031:2050]);
	ba_low_ice_proj(:,ifile)=exp(B1(1) + B1(2)*spi4_3(21:40) + B1(3)*tx_4_9(21:40));
end

%% Low-Ice EC-Earth:
ba_EC_proj=repmat(NaN,20,10);
for ifile=1:10
	load(sprintf('/home/sixto/Dropbox/ivana_ice_fires/review/data/low_ice_ecearth/spi4_3_%d_nostd_baileys.mat',ifile-1));
	load(sprintf('/home/sixto/Dropbox/ivana_ice_fires/review/data/low_ice_ecearth/tx_4_9_%d_nostd_baileys.mat',ifile-1));
    tx_4_9(find(isnan(tx_4_9.*spi4_3)))=NaN;spi4_3(find(isnan(tx_4_9.*spi4_3)))=NaN;
	spi4_3=scale_1981_2000(spi4_3,[1981:2000 2031:2050]);
	tx_4_9=scale_1981_2000(tx_4_9,[1981:2000 2031:2050]);
	ba_EC_proj(:,ifile)=exp(B1(1) + B1(2)*spi4_3(21:40) + B1(3)*tx_4_9(21:40));
end

aux1=nansum(ba_cmip5_proj)';
aux2=nansum(ba_cmip6_proj)';
aux3=nansum(ba_low_ice_proj)';
aux4=nansum(ba_EC_proj)';
tabla=[nanmean(aux1) prctile(aux1,[2.5 50 97.5]);nanmean(aux2) prctile(aux2,[2.5 50 97.5]);nanmean(aux3) prctile(aux3,[2.5 50 97.5]);nanmean(aux4) prctile(aux4,[2.5 50 97.5])];

%% 		       Mean   Prc 2.5    Median  Prc 97.5
%% CMIP5:   33.4485   13.1524   31.8636   65.1695
%% CMIP6:   45.7477   14.6096   40.9300  106.8763
%% Low-Ice:  9.9394    5.9245    9.8815   14.4641
%% EC-Earth: 9.2405    6.6772    9.3004   12.4788
%% 
%% 1951-2000   18.3026
%% 1981-2000   10.7761
%% 2001-2020   28.9918
