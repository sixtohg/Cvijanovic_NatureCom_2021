function dataout=scale_1981_2000(datain,years)
%datain=spi4_3';
%years=1950:2020;
%datain=years;


period=1981:2000;
[~,~,iok]=intersect(period, years);
avg=nanmean(datain(iok),1);
devstd=std(datain(iok),'omitnan');
dataout=(datain-avg)/devstd;


    




     