function [dbz] = dbah(mat,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% CREATED: 2019-11-22
% LAST MODIFIED: 2020-07-22
% db scale with max at zero
% absolute value 
% zero padded hilbert
% reference max optional
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optargin = size(varargin,2);

dbz=hilbert(mat,size(mat,1)*2);
dbz=abs(dbz);
dbz=db(dbz);
if(optargin==0)
    dbz=dbz-maxmax(dbz);
end
if(optargin==1)
dbz=dbz-db(varargin{1});
end

    
    

dim=ndims(mat);
if(dim==1)
    dbz=dbz(1:size(dbz,1)/2);
end
if(dim==2)
    dbz=dbz(1:size(dbz,1)/2,:);
end
if(dim==3)
    dbz=dbz(1:size(dbz,1)/2,:,:);
end
if(dim==4)
    dbz=dbz(1:size(dbz,1)/2,:,:,:);
end
if(dim>=5)
    disp('ERROR dimension too large')
end
