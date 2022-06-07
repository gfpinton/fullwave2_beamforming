function [dbz] = pwah(mat,pwr,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% CREATED: 2022-06-06
% LAST MODIFIED:  2022-06-06
% power scale with max at zero
% absolute value 
% zero padded hilbert
% reference max optional
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optargin = size(varargin,2);

dbz=hilbert(mat,size(mat,1)*2);
dbz=abs(dbz);
dbz=powcompress(dbz,pwr);
if(optargin==0)
    dbz=dbz-maxmax(dbz);
end
if(optargin>=1)
dbz=dbz-powcompress(varargin{1},pwr);
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
