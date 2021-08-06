clear all,clc,close all
workPath=[pwd '/../'];cd([workPath 'script/'])
%%

years=(1951:2020);
alpha = 0.05;
NB=1000;
dir_fire=[workPath 'data/'];
dir_out=[workPath 'paper/figures_temp/'];
%% dir_fire='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/';
%% dir_out='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/paper/figures_temp/';
domain='baileys';
%%
savename = [dir_fire,'fire_1950_2020_',domain,'_mjjas_gt300.mat'];
load(savename,'FIRE')
anni_fire=1950:2020;
[C,IA,IB] = intersect(years,anni_fire);
FIRE=FIRE(IB)/1000;

%% figure; plot(FIRE)
yor=(log(FIRE)');
%% figure; plot(yor)

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

format shortG

%% regression only SPI
Xorig1 = [ones(size(SP,1),1) SP  years' (years.*years)' ]; %osservati
[B1,BINT,R,RINT,STATS1] = regress((yor),Xorig1,alpha);
STATS1
[~,akc]=aic2(length(find(isnan(yor)==0)),size(Xorig1,2),nanvar(R),nanvar((yor)))

%% regression only TSMAX
Xorig1 = [ones(size(SP,1),1) TX  years' (years.*years)' ]; %osservati
[B1,BINT,R,RINT,STATS1] = regress((yor),Xorig1,alpha);
STATS1
[~,akc]=aic2(length(find(isnan(yor)==0)),size(Xorig1,2),nanvar(R),nanvar((yor)))

%% regression all drivers with tsmax time interaction
Xorig1 = [ones(size(SP,1),1) SP TX years' (years.*years)' TX.*years']; %osservati
[B1,BINT,R,RINT,STATS1] = regress((yor),Xorig1,alpha);
STATS1
[~,akc]=aic2(length(find(isnan(yor)==0)),size(Xorig1,2),nanvar(R),nanvar((yor)))

bootdat = [(yor),Xorig1];
bootb1 = bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat);
bootCI1=prctile(bootb1,[2.5 97.5]);

format long
[B1 bootCI1(1,:)'  bootCI1(2,:)']
format long
B1

format shortG
B1


%% regression all drivers and only climate drivers
Xorig1 = [ones(size(SP,1),1) SP TX years' (years.*years)' ]; %osservati
[B1,BINT,R,RINT,STATS1] = regress((yor),Xorig1,alpha);
STATS1
[~,akc]=aic2(length(find(isnan(yor)==0)),size(Xorig1,2),nanvar(R),nanvar((yor)))



format long
B1

format shortG
B1


%% STATS1
%[~,akc]=aic2(length(find(isnan(yor)==0)),size(Xorig1,2),nanvar(R),nanvar(nanpstd(yor)));

Y1 = (Xorig1*B1);
Y2 = B1(1) + B1(2) * SP + B1(3) * TX + B1(4) * (years') + B1(5) * (years.*years)';
Y3 = B1(1) + B1(2) * SP + B1(3) * TX + B1(4) * mean(years') + B1(5) * mean(years.*years)';
Y4 = B1(1) + B1(2) * SP + B1(3) * TX + B1(4) * mean((1981:2000)') + B1(5) * mean((1981:2000).*(1981:2000))';
Y5 = B1(1) + B1(2) * SP + B1(3) * TX + B1(4) ;

figure;
plot(years,(yor),'-ko','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
hold on
%plot(years,(Y1),'k:o','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',6)
%plot(years,(Y2),'r:.','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',6)
plot(years,(Y3),'-go','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',6)
plot(years,(Y4),'-gs','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',6)

mean(exp(yor(31:50)))
mean(exp(Y3(31:50)))
mean(exp(Y4(31:50)))
mean(exp(Y5(31:50)))


[a,b]=corrcoef((yor),Y1);

bootdat = [(yor),Xorig1];
bootb1 = bootstrp(NB, @(x) regress(x(:,1),x(:,2:end)),bootdat);
%% bootb1 = sort(bootb1);
bootCI1=prctile(bootb1,[2.5 97.5]);

format long
[B1 bootCI1(1,:)'  bootCI1(2,:)']
cc
%% bootCI1 = [bootb1(25,:); bootb1(975,:)]  %Report 5% confidence intervals
%% B1
clear bootdat
for ib=1:1000
    all_drivers(:,ib) = bootb1(ib,1) + bootb1(ib,2) * SP + bootb1(ib,3)  * TX + bootb1(ib,4) * years' + bootb1(ib,5) * (years.*years)';
    climate_drivers(:,ib) = bootb1(ib,1) + bootb1(ib,2) * SP + bootb1(ib,3)  * TX + bootb1(ib,4) * mean((1981:2000)')  + bootb1(ib,5) * mean((1981:2000).*(1981:2000))';
end
%%% 
%% hay que averiguar si es mas eficaz en lugar de los que es ahora fig 2d, 
%% una grafica estilo williams et al 2019, es decir, promedio BA/a??o periodo 
%% 1981-2000, promedio BA/a??o periodo 2001-2020, delta (promedio BA/a??o periodo 
%% 2001-2020 - promedio BA/a??o periodo 1981-2000)
%% puedes hacer esta grafica? 
%%% 
auxBar=[mean(exp(yor(31:50))) mean(exp(yor(51:70))) mean(exp(yor(51:70)))-mean(exp(yor(31:50)));...
mean(exp(Y1(31:50))) mean(exp(Y1(51:70))) mean(exp(Y1(51:70)))-mean(exp(Y1(31:50)));...
mean(exp(Y3(31:50))) mean(exp(Y3(51:70))) mean(exp(Y3(51:70)))-mean(exp(Y3(31:50)))];
figure,bar(auxBar','grouped'),hold on,
plot([1 1],prctile(mean(exp(all_drivers(31:50,:)))',[2.5 97.5]),'-k')
plot([2 2],prctile(mean(exp(all_drivers(51:70,:)))',[2.5 97.5]),'-k')
plot([3 3],prctile(mean(exp(all_drivers(51:70,:)))',[2.5 97.5])-prctile(mean(exp(all_drivers(31:50,:)))',[2.5 97.5]),'-k')
plot([1.225 1.225],prctile(mean(exp(climate_drivers(31:50,:)))',[2.5 97.5]),'-k')
plot([2.225 2.225],prctile(mean(exp(climate_drivers(51:70,:)))',[2.5 97.5]),'-k')
plot([3.225 3.225],prctile(mean(exp(climate_drivers(51:70,:)))',[2.5 97.5])-prctile(mean(exp(climate_drivers(31:50,:)))',[2.5 97.5]),'-k')
legend('Obs','All drivers','Climate drivers')
set(gca,'XtickLabel',{'1981-2000','2001-2021','Change'})

figure;
plot(years,(yor),'-ko','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
hold on
% plot(years,(Y1),'k:o','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',6)
plot(years,prctile(climate_drivers,[2.5],2),'k:o','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',6)
plot(years,prctile(climate_drivers,[97.5],2),'k:o','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',6)
legend('Observed','Modelled (lower bound)','Modelled (upper bound)','Location','Best')
xlabel('Years') 

figure;
plot(years,(yor),'-ko','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
hold on
% plot(years,(Y1),'k:o','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',6)
plot(years,prctile(all_drivers,[2.5],2),'k:o','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',6)
plot(years,prctile(all_drivers,[97.5],2),'k:o','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',6)
legend('Observed','Modelled (lower bound)','Modelled (upper bound)','Location','Best')
xlabel('Years') 

figure; plot((yor),Y1,'o')
%% figure
%% figure;drawSpread((dataYM),'xvalues',years','lines','no','indexes',(1:3)','colors',mycol,'widths',widths,'colorsg',colorsg);
indices{1}=[1:size(all_drivers,2)];
boundary=cell(1,1);boundary{1}={'prc2.5','prc97.5'}; 
h=figure;drawSpread(all_drivers,'xvalues',years','indexesg',indices,'colorsg',[255 0 0]/255,'lines','no','boundary',boundary,'alphasg',[0.5 0.5]); hold on
drawSpread(climate_drivers,'xvalues',years','indexesg',indices,'colorsg',[255,165,0]/255,'lines','no','boundary',boundary,'alphasg',[0.5 0.5]);


indices{1}=[1:size(all_drivers,2)];indices{2}=[1:size(all_drivers,2)];
boundary=cell(2,1);boundary{1}={'prc2.5','prc97.5'}; boundary{2}={'prc2.5','prc97.5'};
h=figure;drawSpread(all_drivers,'xvalues',years','indexesg',indices,'colorsg',[255 0 0;255 0 0]/255,'lines','no','boundary',boundary,'alphasg',[0.5 0.5]); hold on
drawSpread(climate_drivers,'xvalues',years','indexesg',indices,'colorsg',[255,165,0;255,165,0]/255,'lines','no','boundary',boundary,'alphasg',[0.5 0.5]);
%plot(years,(Y1'),'k:o','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',6)
plot(years,(yor'),'k-o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
legend('All drivers (2.5-97.5)','Climate drivers (2.5-97.5)','','','Observed','Location','Best')
xlabel('Years','FontSize',18) 
ylabel('log(BA) [log(1000 km^2)]','FontSize',18) 
%gridxy(years(ith),'Color','k','Linestyle',':') ;
set(gca,'FontSize',18,'xlim',[years(1)-1 years(end)+1])
set(gcf, 'PaperPositionMode', 'auto','renderer', 'painters');
file=[dir_out,'figure1c_temp_drivers.eps']
print( gcf, '-depsc2', file ,'-painters')

set(gcf,'PaperType','A4')

print( gcf, '-dpdf', [file(1:end-4),'.pdf'] ,'-painters')

% 



%% as in calmanti et al. 2007

%years=years';
alpha = 0.05;

oos2=[]; 
  
p=1;
numanni=length(years);
Xorig = [ones(size(SP,1),1) SP TX years' (years.*years)' ]; %osservati
for i=1:p:numanni;
    itest=[i:i+(p-1)];
    itrain1=setdiff([1:numanni],itest);
    %X = [ones(size(PSM1(itrain1))) PSM1(itrain1) PSM2(itrain1) TMAX1(itrain1) TMIN1(itrain1) TMIN2(itrain1)];
    X=Xorig(itrain1,:); 
    
    y=yor(itrain1);
    [B,BINT,R,RINT,STATS] = regress(y,X,alpha);
    sigma=std(R);
    %X = [ones(size(PSM1(itest))) PSM1(itest) PSM2(itest) TMAX1(itest) TMIN1(itest) TMIN2(itest)];
    clear X
    X=Xorig(itest,:); 
   
    %plot(years,nf,'r-')
    %hold on
    oos1=[]; 
    for iboot=1:1000
        noise=normrnd(0,sigma,length(itest),1);
        Y = (X*B)+noise;
        oos1=[oos1,Y];
        %plot(years(itest),exp(Y)','k:')
    end   
    oos2=[oos2;oos1];
    %plot(years,nf,'r-') 
    Boos(i,:)=B;
    Boos(i,:)=B;
end    

for i=1:numanni
    %b2(i,:)=[min(exp(oos2(i,:))-nf(i)),max(exp(oos2(i,:))-nf(i))];
    Ymed(i)=prctile((oos2(i,:)),50);
    b1(i,:)=[prctile((oos2(i,:)),2.5),prctile((oos2(i,:)),97.5)];
    
end 

figure1 = figure('PaperSize',[20.98 29.68]);
%ylim([-4 4])
%xlim([-4 4])
% Create axes
axes1 = axes('Parent',figure1);
hold('all');

neg=b1(:,1);
pos=b1(:,2);

for ip=1:length(Ymed)
    ly=neg(ip);
    uy=pos(ip);
    %l1=line([log(BAS(ip)) log(BAS(ip))],[ly uy],'color',[128 128 128]/255,'LineWidth',0.5);
    l1=line([yor(ip) yor(ip)],[ly uy],'color',[128 128 128]/255,'LineWidth',0.5);
    
    %l2=line([TMEAN_ECO_YM(ip)-0.1 TMEAN_ECO_YM(ip)+0.1],[ly ly],'color',[128 128 128]/255,'LineWidth',0.5);
    %l3=line([TMEAN_ECO_YM(ip)-0.1 TMEAN_ECO_YM(ip)+0.1],[uy uy],'color',[128 128 128]/255,'LineWidth',0.5);
end
%scatter1 =scatter(log(BAS),Ymed,150,years,'MarkerFaceColor','flat','MarkerEdgeColor',[0 0 0])
scatter1 =scatter(yor,Ymed,150,years,'MarkerFaceColor','flat','MarkerEdgeColor',[0 0 0])
h=colorbar
colormap('autumn')
scatter1.MarkerFaceAlpha = 0.8;
title(h,'Year','FontSize',18)
%line([-4 4],[-4 4],'color','k')
line([-5 3.5],[-5 3.5],'color','k')


xlim([-5 3.5])
ylim([-5 3.5])
xlabel('Observed log(BA) [log(1000 km^2)]','FontSize',12)
ylabel('Predicted log(BA) [log(1000 km^2)]','FontSize',12)
set(gca,'FontSize',12)

set(gcf,'PaperType','A4')

file=[dir_out,'figure1d_temp.pdf']
print( gcf, '-dpdf', file ,'-painters')



file=[dir_out,'figure1d_temp.eps']
print( gcf, '-depsc2', file ,'-painters')



%% Create figure
st=1;
for i=1:size(years,1)
    Ymed(i)=prctile((oos2(i,:)),50);
    b1(i,:)=[prctile((oos2(i,:)),2.5)-Ymed(i),prctile((oos2(i,:)),97.5)-Ymed(i)];
    b2(i,:)=[prctile((oos2(i,:)),25)-Ymed(i),prctile((oos2(i,:)),75)-Ymed(i)];
end    
b1=abs(b1);
b2=abs(b2);

[a,b]=corrcoef((yor),Ymed,'rows','complete')
 



dataYM(:,1)=Ymed';
dataYM(:,2)=(Ymed'-b1(:,1));
dataYM(:,3)=(Ymed'+b1(:,2));

dataYM2(:,1)=Ymed';
dataYM2(:,2)=(Ymed'-b2(:,1));
dataYM2(:,3)=(Ymed'+b2(:,2));
colorsg    = [1, 0, 0]; %red
colorsg2    = [255 69 0]/255; %orange
mycol= [0 0 0;1 0 0; 1 0 0]; %nero, red red
mycol2= [0 0 0;[255 69 0]/255; [255 69 0]/255]; %nero, orange orange red


%colorsg    = [0.85, 0.85, 0.85]; %grey
%mycol= [0 0 0;0.4 0.4 0.4; 0.4 0.4 0.4]; %nero, grigio grigio
widths= [2 1 1]';
numanni=length(years); %totale lunghezza serie
ith=1:st:numanni;


%% figure;drawSpread((dataYM),'xvalues',years','lines','no','indexes',(1:3)','colors',mycol,'widths',widths,'colorsg',colorsg);
indices{1}=[1:size(oos2,2)];indices{2}=[1:size(oos2,2)];
boundary=cell(2,1);boundary{1}={'prc2.5','prc97.5'};boundary{2}={'prc25','prc75'};
h=figure;drawSpread(oos2,'xvalues',years','indexesg',indices,'colorsg',[255 0 0;101 101 101]/255,'lines','no','boundary',boundary,'alphasg',[0.5 0.5]);hold on,
plot(years,(yor'),'k-o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
plot(years,(Ymed'),'k:o','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',6)
legend('Modelled (2.5-97.5)','Modelled (25-75)','Observed','Location','Best')
xlabel('Years','FontSize',18)
ylabel('log(BA)','FontSize',18)
%gridxy(years(ith),'Color','k','Linestyle',':') ;
set(gca,'FontSize',18,'xlim',[years(1)-1 years(end)+1])
set(gcf, 'PaperPositionMode', 'auto','renderer', 'painters');
file=[dir_out,'figure1c_temp.eps']
print( gcf, '-depsc2', file ,'-painters')


%% hold on
%% plot(years,(yor'),'k-o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
%% plot(years,(Ymed'),'k:o','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',6)
%% legend('Modelled','Observed','Location','Best')
%% xlabel('Years','FontSize',18) 
%% ylabel('log(BA)','FontSize',18) 
%% %gridxy(years(ith),'Color','k','Linestyle',':') ;
%% set(gca,'FontSize',18) 
%% file=[dir_out,'OS_',num2str(st),'_nfsu.eps']
%% set(gcf, 'PaperPositionMode', 'auto','renderer', 'painters');




% figure;drawSpread((dataYM),'xvalues',years','lines','no','indexes',(1:3)','colors',mycol,'widths',widths,'colorsg',colorsg);
% hold on
% plot(years,(yor'),'k-o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6)
% plot(years,(Ymed'),'k:o','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',6)
% %ylim([-3.5 3])
% legend('Modelled','Observed','Location','NorthWest')
% %drawSpread((dataYM2),'xvalues',years','lines','no','indexes',(1:3)','colors',mycol2,'widths',widths,'colorsg',colorsg2);
% 
% xlabel('Years','FontSize',18) 
% ylabel('log(BA) [log(1000 km^2)] ','FontSize',18) 
% %gridxy(years(ith),'Color','k','Linestyle',':') ;
% set(gca,'FontSize',18) 
% file=[dir_out,'OS_',num2str(st),'_nfsu.eps']
% set(gcf, 'PaperPositionMode', 'auto','renderer', 'painters');
% xlim([1949 2022])
% file=[dir_out,'figure1c_temp.eps']
% print( gcf, '-depsc2', file ,'-painters')
