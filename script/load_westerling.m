clear all
close all
%%

years=1980:2004;
%%
dir_fire='/Users/marco/Documents/virtualBox/lavori/drought_fire_cal/';
filename = [dir_fire,'ALL.ACRE.8004.xlsx'];
[NUM,TXT,RAW]=xlsread(filename);
clear TXT RAW

coord(:,1)=NUM(2,3:end)';
coord(:,2)=NUM(1,3:end)';

drawStations(coord,'resolution','high')


%% fire data
% 1000 acres = 4.04685642 square kilometers
% 1 acre = 4.04685642/1000 square kilometers
BAm_all=NUM(3:end,3:end)*(4.04685642/1000); %acres to kmq
%figure;plot(BAm)


%--> it seems that there are also fires < 1000 acres
% thus I exclude NUM(:,1)<=4.04685642)
%ino=find(NUM<=4.04685642);
%NUM(ino,:)=[];



%% load california
filename = '/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/shapefiles/baileys/Baileys_ecoregions_cal_wgs84.shp';
S = shaperead(filename);
geoshow(S)

for izone=1:size(S,1)
    dum=S(izone).X';
    xv=[dum;dum(1)]; clear dum
    dum=S(izone).Y';
    yv=[dum;dum(1)]; clear dum
    in = inpolygon(coord(:,1),coord(:,2),xv,yv);
    figure;drawStations(coord(in,:))

    %% monthly aggregation
    
    k=0;
    for iyear=years
        for im=1:12
            k=k+1;
            BAm(k,izone)=nansum(BAm_all(k,in),2);
        end
    end
end



for izone=1:size(S,1)
    dum=nansum(BAm(:,izone),2);
    for im=1:12
        BAmm(im,:)=dum(im:12:end,:);
    end
    
    color_fire=[255,69,0]/255;
    figure;
    aboxplot2(BAmm','colormap',color_fire,'labels',[1:12]); % Advanced box plot
    xlabel('Months','FontSize',20)
    ylabel('Burned Area (km^2)','FontSize',20)
    set(gca,'FontSize',18)
    title('BA','FontSize',24)
    %ylim([0 4100])
    %file=[dir_out,'annual_cycle_BA_firep14_2.eps']
    %print( gcf, '-depsc2', file ,'-painters')
end

%% JJASO aggregation
BAs=zeros(length(years),size(S,1));
for izone=1:size(S,1)
    
    for iyear=1:length(years)
        i1=(iyear-1)*12+5
        i2=i1+4
        BAs(iyear,izone)=nansum(BAm(i1:i2,izone));
        
    end
    
    
    
    x1=years;
    y1=BAs(:,izone);
    figure;
    hl1 = bar(x1,y1,'k');
    xlim([x1(1)-1 x1(end)+1]);
    %legend('BA','Location','NorthWest')
    ax1 = gca;
    set(ax1, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     )
    xlabel('Years')
    set(get(ax1,'Ylabel'),'String','Burned Area (km^2)','color','k')
    ylimits1 = get(ax1,'YLim');
    %set(ax1,'YTick',[ylimits1(1):1000:ylimits1(2)])
    %file=[dir_out,'mtbs_AB_JJA_firep14_2.eps']
    %print( gcf, '-depsc2', file ,'-painters')
    
    %[pValue,trend]=testTrend_MT_sen(BAs)
end



%% save JJASO fire
dir_out='/Users/marco/Documents/dati/fire_us/';
savefile = [dir_out,'westerling_1980-2004_mjjas.mat'];
save(savefile,'BAs')


