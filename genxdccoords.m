function [incoords outcoords incoords2 outcoords2]=genxdccoords(nE,ptch,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: 2021-02-04
% LAST MODIFIED: 2021-04-13
% Generate transducer coordinates with element labels
% nE: number of elements
% ptch: pitch in pixels
% nX: centering variable to center around nX/2
% kerf: kerf width in pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
optargin=size(varargin,2);

kerf=0;
if(optargin==2)
  kerf=varargin{2}
end

nEpx=(nE*(ptch-kerf)+(nE-1)*kerf); % total pixel width of array
incoords=zeros(nEpx*3,3);
for k=1:3
  for i=1:nE
    for j=1:ptch-kerf
      incoords((i-1)*ptch+j+(k-1)*nEpx,1)=(i-1)*ptch+j; % 
      incoords((i-1)*ptch+j+(k-1)*nEpx,2)=k; % thickness of array 
      incoords((i-1)*ptch+j+(k-1)*nEpx,3)=i; % element number
    end
    if(i<nE)
      for j=ptch-kerf+1:ptch
	incoords((i-1)*ptch+j+(k-1)*nEpx,1)=(i-1)*ptch+j; % 
	incoords((i-1)*ptch+j+(k-1)*nEpx,2)=k; % thickness of array 
	incoords((i-1)*ptch+j+(k-1)*nEpx,3)=0; % element number
      end
    end
  end
end

outcoords=zeros((nE*(ptch-kerf)+(nE-1)*kerf),3);
k=4
for i=1:nE
  for j=1:ptch-kerf
    outcoords((i-1)*ptch+j,1)=(i-1)*ptch+j; % 
    outcoords((i-1)*ptch+j,2)=k; % thickness of array 
    outcoords((i-1)*ptch+j,3)=i; % element number
  end
  if(i<nE)
    for j=ptch-kerf+1:ptch
      outcoords((i-1)*ptch+j,1)=(i-1)*ptch+j; % 
      outcoords((i-1)*ptch+j,2)=k; % thickness of array 
      outcoords((i-1)*ptch+j,3)=0; % element number
    end
  end
end

for tt=1:nE
  idt=find(outcoords(:,3)==tt);
  outcoords2(tt,:)=mean(outcoords(idt,:),1);
  idt=find(incoords(:,3)==tt);
  incoords2(tt,:)=mean(incoords(idt,:),1);
end

if(optargin>=1) % center around nX/2
  nX=varargin{1};
  incoords(:,1)=round(incoords(:,1)-mean(incoords(:,1))+nX/2);
  outcoords(:,1)=round(outcoords(:,1)-mean(outcoords(:,1))+nX/2);
  incoords2(:,1)=round(incoords2(:,1)-mean(incoords2(:,1))+nX/2);
  outcoords2(:,1)=round(outcoords2(:,1)-mean(outcoords2(:,1))+nX/2);
end
  

  
