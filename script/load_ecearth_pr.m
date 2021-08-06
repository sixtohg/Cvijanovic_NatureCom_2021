function [datagrid]=load_ecearth_pr(filename)


% Example: 
%

%% Default parameters

[dir_dro, baseFileName, ~] = fileparts(filename);
savefile = sprintf('%s.mat', baseFileName);

sf=1; %scale factor
varname='pr';
savefile=[dir_dro '/' savefile];

%% nc to mat

ncid = netcdf.open(filename,'NC_NOWRITE');
[numdims,nvars,~] = netcdf.inq(ncid);


% for ii = 0:natts-1
%     fieldname = netcdf.inqAttName(ncid, netcdf.getConstant('NC_GLOBAL'), ii);
%     fileinfo.(fieldname) = netcdf.getAtt(ncid,netcdf.getConstant('NC_GLOBAL'), fieldname );
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
    datagrid = netcdf.getVar(ncid,ii-1);
    
    if (isfield(tmpstruct, 'FillValue') )
        datagrid( datagrid == tmpstruct.FillValue ) = NaN;
    end
    %if( isfield(tmpstruct, 'scale_factor') )
    %    tmpstruct.scale_factor
    %    data = double(data) * tmpstruct.scale_factor;
    %end
    %if( isfield(tmpstruct, 'add_offset') )
    %    data = data + tmpstruct.add_offset;
    %end
    if( isnumeric(datagrid) && ndims(datagrid) > 2 )
        datagrid = permute(datagrid, [2 1 3:ndims(datagrid)]);
    elseif ( isnumeric(datagrid) && ndims(datagrid) == 2 )
        datagrid = datagrid';
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
    
    if strcmp(name,varname)
        dum=datagrid;
    end
    
end
netcdf.close(ncid);



size(dum)
dum2=dum;
dum2(:,(1:size(dum,2)/2),:)=dum(:,(size(dum,2)/2+1:size(dum,2)),:);
dum2(:,(size(dum,2)/2+1:size(dum,2)),:)=dum(:,(1:size(dum,2)/2),:);

%# dum[1:(nrow(obs)/2),,]=obs[(nrow(obs)/2+1):nrow(obs),,]
%  # dum[(nrow(obs)/2+1):nrow(obs),,]=obs[1:(nrow(obs)/2),,]

tmp = NaN*zeros(size(dum,3),size(dum,1)*size(dum,2));

%Passo a matrix(days,space)
for im=1:size(dum,3)
    k=0;
    for j=1:size(dum,2)
        for i=1:size(dum,1)
            
            k=k+1;
            tmp(im,k)=dum2(i,j,im);
        end
    end
end
clear dum

  
datagrid=tmp*sf;
clear tmp dum
%% salva i file

%data(data<-999)=NaN;
save(savefile,'datagrid')

%save(filename,'dumvar','-v7.3');

coord=load_coord('ecearth');
pr=nanmean(datagrid,1);
drawStations(coord,'size',1,'resolution','high','israster','true','color',pr','colormap',(jet))
%caxis([-3 3])
h=colorbar


