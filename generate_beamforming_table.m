function [idps_mat] = generate_beamforming_table(lats,deps,xducercoords,nT,dY,dT,c0,fnumber)

  idps_mat=zeros(length(lats),length(deps),nT+1);
   
  for ii=1:length(lats)
    lat=lats(ii);
    for jj=1:length(deps)
      dep=deps(jj);
      fcen=([lat/dY dep/dY]);
      idx=find(abs(xducercoords(:,1)-fcen(1))<=fcen(2)/fnumber);
      dd=focusProfile(fcen,xducercoords(idx,:),dT/dY*c0);
      idt=2*dep/double(c0)/(dT); %% NO IDT0
 
      idp=double((nT*(idx-1))+idt+dd); %no rounding
	
      idps_mat(ii,jj,1)=length(idp);
      idps_mat(ii,jj,2:length(idp)+1)=(idp);
    end
  end
  

  
