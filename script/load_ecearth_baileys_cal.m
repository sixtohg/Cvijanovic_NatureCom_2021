clear all
close all

%% Default parameters
% DIR='mac'
%[dir_oss,dir_out]=choose_dir(DIR);
dir_oss='/Users/marco/Documents/dati/ivana_ice_fires/Sixto_ECEarth3/low_ice/';

%% files
list_files_pr = dir(fullfile(dir_oss, '*pr_allmon_y1030.nc'));
list_files_tasmax = dir(fullfile(dir_oss, '*tasmax_allmon_y1030.nc'));

%% ecoregions
filename_pr = '/Users/marco/Dropbox/estcena/scripts/ivana_ice_fires/review/data/shapefiles/baileys/Baileys_ecoregions_cal_wgs84.shp';
eco = shaperead(filename_pr, 'UseGeoCoords', true); % imports lat/lon data

%geoshow(eco)

%% ecearth coord

coord=load_coord('ecearth');

%% load tasmax and pr


for ircm=1:length(list_files_pr)
    list_files_pr(ircm).name
    list_files_tasmax(ircm).name
    
    
    filename_pr=[dir_oss list_files_pr(ircm).name];
    [~, baseFileName_pr, ~] = fileparts(filename_pr);
    savefile_pr = sprintf('%s.mat', baseFileName_pr);
    pr=load_ecearth_pr(filename_pr);
    
    filename_tasmax=[dir_oss list_files_tasmax(ircm).name];
    [~, baseFileName_tasmax, ~] = fileparts(filename_tasmax);
    savefile_tasmax = sprintf('%s.mat', baseFileName_tasmax);
    tasmax=load_ecearth_tasmax(filename_tasmax);
    
    dum=eco.Lon';
    xv=[dum;dum(1)]; clear dum
    dum=eco.Lat';
    yv=[dum;dum(1)]; clear dum
    
    in = inpolygon(coord(:,1),coord(:,2),xv,yv);
    
    figure;drawStations(coord(in,:));
    
    pr_eco=nanmean(pr(:,in),2)*86400; %kg m-2 s-1 --> mm/day
    tasmax_eco=nanmean(tasmax(:,in),2)-273.15; %kelvin --> celsius
    
    savename=[dir_oss '' savefile_pr]
    save(savename,'pr_eco')
    
    savename=[dir_oss '' savefile_tasmax]
    save(savename,'tasmax_eco')
    
    clear pr_eco tasmax_eco
    
    close all
    
end