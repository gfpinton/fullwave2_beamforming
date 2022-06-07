function [focs bmws] = genlinearfocs(nTx,bmw,dep,apex)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: 2021-02-04
% LAST MODIFIED: 2022-05-10
% Generate centered focal delays for linear scan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:nTx
  bmws(n)=bmw*(n-(nTx+1)/2);
  focs(n,1)=bmw*(n-(nTx+1)/2);
  focs(n,2)=dep;
end
focs(:,1)=focs(:,1)+apex(1);
focs(:,2)=focs(:,2)+apex(2);
