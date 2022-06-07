function [idps_mat_rt] = generate_beamforming_receive_table(lats,deps,xducercoords,nT,dY,dT,c0,fnumber)

  idps_mat_rt=zeros(length(lats),length(deps),nT+1); % check this size
   
  for ii=1:length(lats)
    lat=lats(ii);
    for jj=1:length(deps)
      dep=deps(jj);
      fcen=([lat/dY dep/dY]);
      idx=find(abs(xducercoords(:,1)-fcen(1))<=fcen(2)/fnumber);
      [dd mdd]=focusProfile2(fcen,xducercoords(idx,:),dT/dY*c0);
      idp=double((nT*(idx-1))+mdd); %no rounding
	
      idps_mat_rt(ii,jj,1)=length(idp);
      idps_mat_rt(ii,jj,2:length(idp)+1)=(idp);
    end
  end
  
