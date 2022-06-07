function [icmat] = focusCoordsDD (dd,icvec,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: JAN 4, 2021
% Focus coordinates
% mdd: offset in time pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  mdd=min(min(dd));
  optargin=size(varargin,2);
  if(optargin==1)
    if(varargin{1}<mdd)
      mdd=varargin{1};
    else
      disp(['Error in focusCoordsDD, setting mdd to ' num2str(mdd)])
  end
end
dd = dd-mdd;

icmat = zeros(length(dd),length(icvec));
for i=1:length(dd)
  icmat(i,dd(i)+1:end)  = icvec(1:end-dd(i));
end
