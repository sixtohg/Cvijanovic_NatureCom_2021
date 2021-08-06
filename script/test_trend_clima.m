
%spi not significan decreasing trend

clear all
close all

%%

%years=(1981:2020);
years=1951:2020;
alpha = 0.05;
NB=1000;
dir_fire='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/';
dir_out='/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/paper/figures_temp/';
domain='baileys';

%%
filename = [dir_fire,'spi4_3_nclimgrid.mat'];
load(filename);
filename = [dir_fire,'tx_4_9_nclimgrid.mat'];
load(filename);


%% linear trend

x=years;
y=spi4_3';

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

[a,b]=testTrend_MT_sen(y')

%% quadratic trend



n=2;
p = polyfit(x,y,n)
f = polyval(p,x);

dum=corrcoef(y,f);
R2=dum(1,2)^2;
R=y-f;
[ak,akc]=aic2(length(y),size(p,2),var(R),var(y));
 
 
figure; plot(years,y,'o',years,f,'-')
legend('log(NF)',['y=' num2str(p(1), '%+10.4f') 'x^2'  ...
                       num2str(p(2), '%+10.4f') 'x'  ...
                       num2str(p(3), '%+10.4f')   ...
       '; R^2=' num2str(R2, '%10.2f')  '; AIC_c=' num2str(akc, '%10.2f') ],'Location','Best')
xlabel('Years'); 
ylabel('log(BA)');



%% cubic trend



n=3;
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




