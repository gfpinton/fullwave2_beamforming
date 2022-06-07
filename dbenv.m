function [dbz] = dbenv(mat,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% CREATED: 2020-07-26
% LAST MODIFIED: 2020-07-26
% db scale with max at zero
% absolute value 
% zero padded envelope
% reference max optional
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optargin = size(varargin,2);

envtap=60;
if(optargin>=2)
    envtap=varargin{2};
end
mat2=mat; mat2(end+1:2*end,:)=0;

dbz=envelope(mat2,envtap);
dbz=db(dbz);
if(optargin==0)
    dbz=dbz-maxmax(dbz);
end
if(optargin>=1)
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
