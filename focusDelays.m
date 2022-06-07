function [dds] = focusDelays (focs,coords,cfl,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: NOV 13, 2013
% LAST MODIFIED: 2022-04-22
% Focus coordinates
% mdd: offset in time pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if(size(focs,2)==2)
    disp('2D')
    dds=zeros(size(coords,1),size(focs,1));
    for n=1:size(focs,1)
      idy=focs(n,1); idz=focs(n,2);
      dd = sqrt((coords(:,1)-idy).^2+(coords(:,2)-idz).^2);
      dd = -round(dd/cfl);
      mdd=min(dd);
      
      optargin=size(varargin,2);
      if(optargin==1)
	mdd=varargin{1};
      end
      dd = dd-mdd;
      dds(:,n)=dd;
    end
  end
  if(size(focs,2)==3)
    disp('3D')
    dds=zeros(size(coords,1),size(focs,1));
    for n=1:size(focs,1)
      idi=focs(n,1); idy=focs(n,2); idz=focs(n,3);
      dd = sqrt((coords(:,1)-idi).^2+(coords(:,2)-idy).^2+(coords(:,3)-idz).^2);
      dd = -round(dd/cfl);
      mdd=min(dd);
      
      optargin=size(varargin,2);
      if(optargin==1)
	mdd=varargin{1};
      end
      dd = dd-mdd;
      dds(:,n)=dd;
    end
  end
