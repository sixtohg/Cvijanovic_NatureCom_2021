clear all
close all

%%

%years=(1981:2020);
years=1950:2020;
alpha = 0.05;
NB=1000;
dir_fire='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/';
dir_out='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/paper/figures_temp/';
domain='baileys';

%%
%savename = [dir_fire,'fire_1950_2020_',domain,'_mjjas.mat'];
savename = [dir_fire,'fire_1950_2020_',domain,'_mjjas_gt300.mat'];
load(savename,'FIRE')
%anni_fire=1950:2020;
%[C,IA,IB] = intersect(years,anni_fire);
%FIRE=FIRE(IB);

FIRE(1)=[];
years(1)=[];

figure; plot(FIRE)


yor=(log(FIRE));

figure; plot(yor)

%% linear trend

x=years;
y=log(FIRE);

n=1;
p = polyfit(x,y,n)
f = polyval(p,x);

dum=corrcoef(y,f);
R2=dum(1,2)^2;
R=y-f;
[ak,akc]=aic2(length(y),size(p,2),var(R),var(y));
 
 
figure; plot(years,y,'o',years,f,'-')
legend('log(BA)',['y=' num2str(p(1), '%+10.4f') 'x' num2str(round(p(2)), '%+10.1f') ...
       '; R^2=' num2str(R2, '%10.2f')  '; AIC_c=' num2str(akc, '%10.2f') ],'Location','Best')
xlabel('Years'); 
ylabel('log(BA)');

%% quadratic trend

x=years;
y=log(FIRE);

n=2;
p = polyfit(x,y,n)
f = polyval(p,x);

dum=corrcoef(y,f);
R2=dum(1,2)^2;
R=y-f;
[ak,akc]=aic2(length(y),size(p,2),var(R),var(y));
 
 
figure; plot(years,y,'o',years,f,'-')
legend('log(BA)',['y=' num2str(p(1), '%+10.4f') 'x^2'  ...
                       num2str(p(2), '%+10.4f') 'x'  ...
                       num2str(p(3), '%+10.4f')   ...
       '; R^2=' num2str(R2, '%10.2f')  '; AIC_c=' num2str(akc, '%10.2f') ],'Location','Best')
xlabel('Years'); 
ylabel('log(BA)');



%% cubic trend

x=years;
y=log(FIRE);

n=3;
p = polyfit(x,y,n)
f = polyval(p,x);

dum=corrcoef(y,f);
R2=dum(1,2)^2;
R=y-f;
[ak,akc]=aic2(length(y),size(p,2),var(R),var(y));
 
 
figure; plot(years,y,'o',years,f,'-')
legend('log(BA)',['y=' num2str(p(1), '%+10.4f') 'x^3'  ...
                       num2str(p(2), '%+10.4f') 'x^'  ...
                       num2str(p(3), '%+10.4f') 'x'  ...
                       num2str(p(4), '%+10.4f')   ...
       '; R^2=' num2str(R2, '%10.2f')  '; AIC_c=' num2str(akc, '%10.2f') ],'Location','Best')

   xlabel('Years'); 
ylabel('log(BA)');


cc

%% califonia nclimgrid
filename = [dir_fire,'nclimgrid_baileys_cal.mat'];
load(filename);
%nclimgrid spi1,2,3,4,5,6,9,12,24,tasmax.tasmin,tas
spi6_3 = nclimgrid(6,3:12:end);
spi6_3(1:2)=[]; %delete years 1948 and 1949

spi3_2 = nclimgrid(3,2:12:end);
spi3_2(1:2)=[]; %delete years 1948 and 1949

spi4_3 = nclimgrid(4,3:12:end);
spi4_3(1:2)=[]; %delete years 1948 and 1949
spi4_3(1)=[]; %delete year 1950

%% PREC vs BA
SP=spi4_3;

figure1 = figure('PaperSize',[20.98 29.68]);
axes1 = axes('Parent',figure1);
hold('all');
scatter1 =scatter(detrend(SP(1:30)),detrend(yor(1:30)),150,years(1:30),'MarkerFaceColor','flat','MarkerEdgeColor',[0 0 0])
colorbar
colormap('autumn')
axesLimits1 = xlim(axes1);
xplot1 = linspace(axesLimits1(1), axesLimits1(2));
fitResults1 = polyfit(detrend(SP(1:30)),detrend(yor(1:30)), 1);
yplot1 = polyval(fitResults1, xplot1);
fitLine1 = plot(xplot1,yplot1,'k:');
ylabel('detrended log(BA)','FontSize',20);
xlabel('detrended $SPI6_{3}$','FontSize',20);
set(gca,'FontSize',18)
%file=[dir_out,'scatter_ba_spi6_3_1950_1980.eps']
%print( gcf, '-depsc2', file ,'-painters')
%saveas(gcf,[file(1:end-4),'.fig'])


figure1 = figure('PaperSize',[20.98 29.68]);
axes1 = axes('Parent',figure1);
hold('all');
scatter1 =scatter(detrend(SP(31:end)),detrend(yor(31:end)),150,years(31:end),'MarkerFaceColor','flat','MarkerEdgeColor',[0 0 0])
colorbar
colormap('autumn')
axesLimits1 = xlim(axes1);
xplot1 = linspace(axesLimits1(1), axesLimits1(2));
fitResults1 = polyfit(detrend(SP(31:end)),detrend(yor(31:end)), 1);
yplot1 = polyval(fitResults1, xplot1);
fitLine1 = plot(xplot1,yplot1,'k:');
ylabel('detrended log(BA)','FontSize',20);
xlabel('detrended $SPI6_{3}$','FontSize',20);
set(gca,'FontSize',18)
%file=[dir_out,'scatter_ba_spi6_3_1981_2018.eps']
%print( gcf, '-depsc2', file ,'-painters')
%saveas(gcf,[file(1:end-4),'.fig'])
years(1:30)
years(31:end)
[rho,pval]=corr(detrend(SP(1:30))',detrend(yor(1:30))')
[rho,pval]=corr(detrend(SP(31:end))',detrend(yor(31:end))')
