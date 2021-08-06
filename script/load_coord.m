function [coord]=load_coord(dataset,varargin)

% Estrae le coordinate presenti in un file netcdf
%
%  Input arguments ([]'s are optional): 
%  dataset     : name of the dataset (e.g.'eobs')    
%
%  Output arguments:
%  coord: n points x 2 (1' column: longitude; 2' column latitude)
%
%  varargin    : optional parameters
%      'DIR'    -   %  home or univ (mac by default)
%
% Example: 
% [coord]=load_coord('ictp','DIR','mac');
%

%% Default parameters

DIR   = 'mac';  % macchina su cui lavoro


%% Read parameters
i=1;
while i<=length(varargin), 
  argok = 1;
  switch varargin{i},
     case 'DIR',      i=i+1; DIR = varargin{i}; 
  otherwise
     argok=0;
  end
  if ~argok, 
    disp(['Ignoring invalid argument #' num2str(i+1)]); 
  end
  i = i+1; 
end


%% scelta dei nomi
dir_dro='/Users/marco/Documents/dati/ivana_ice_fires/Sixto_ECEarth3/';


switch lower(dataset)
    case {'ecearth'};
        nomelon='longitude';
        nomelat='latitude';
        filename=[dir_dro 'r12i1p2f1.20570101_untar.pr_allmon_y1030.nc']
        savefile=[dir_dro 'coord_ecearth.mat']
    otherwise
        disp('Unknown dataset')
end


%% nc to mat

ncid = netcdf.open(filename,'NC_NOWRITE');
[numdims,nvars,natts] = netcdf.inq(ncid);


% for ii = 0:natts-1
%     ii
%     fieldname = netcdf.inqAttName(ncid, netcdf.getConstant('NC_GLOBAL'), ii)
%     fileinfo.(fieldname) = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'), fieldname )
% end

dimension = repmat(struct('name', '', 'length', 0), numdims, 1);

for ii = 1:numdims
    [dimension(ii).name, dimension(ii).length] = netcdf.inqDim(ncid,ii-1);
    
    padlength = min(0, length(dimension(ii).name));
    name_padded = [dimension(ii).name repmat(' ', padlength+1)];
    
    fprintf('%s\t\t%d\n', name_padded, dimension(ii).length)
end



for ii = 1:nvars
    [name, ~, ~, natts] = netcdf.inqVar(ncid,ii-1);
    tmpstruct = struct();
    for jj = 1:natts
        fieldname = netcdf.inqAttName(ncid, ii-1, jj-1);
        if strcmp(fieldname,'_FillValue')
            fieldname2 ='FillValue';
            tmpstruct.(fieldname2) = netcdf.getAtt(ncid, ii-1, fieldname );
        end
        if strcmp(fieldname,'_CoordinateAxisType')
            fieldname2 ='CoordinateAxisType';
            tmpstruct.(fieldname2) = netcdf.getAtt(ncid, ii-1, fieldname );
        end
    end
    data = netcdf.getVar(ncid,ii-1);
    if (isfield(tmpstruct, 'missing_value') )
        data( data == tmpstruct.missing_value ) = NaN;
    end
    if( isfield(tmpstruct, 'scale_factor') )
        data = double(data) * tmpstruct.scale_factor;
    end
    if( isfield(tmpstruct, 'add_offset') )
        data = data + tmpstruct.add_offset;
    end
     if( isnumeric(data) && ndims(data) > 2 )
     data = permute(data, [2 1 3:ndims(data)]);
     elseif ( isnumeric(data) && ndims(data) == 2 )
     data = data';
     end
    if(strcmp('W-E', name))
        name='WE'
    end
    if(strcmp('S-N', name))
        name='SN'
    end
    
    %varinfoname = [name '_info'];
    %assignin('caller', varinfoname, tmpstruct);
    %assignin('caller', name, data);
    
    if strcmp(name,nomelon)
            Longitude=data;
    end
     
    if strcmp(name,nomelat)
            Latitude=data;
    end
        
    name
end

size(Latitude)
size(Longitude)
%Latitude=Latitude';
%Longitude=Longitude';

coord = NaN*zeros(size(Latitude,1)*size(Latitude,2),2);

%Passo a matrix(days,space)

k=1;
for i=1:size(Latitude,2)
    for j=1:size(Latitude,1)
        coord(k,1)=Longitude(j,i)-180;
        coord(k,2)=Latitude(j,i);
        k=k+1;
    end
end





figure;drawStations(coord)

netcdf.close(ncid);


%% salva i file
save(savefile,'coord')    
