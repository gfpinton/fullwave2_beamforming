function [idps_mat_ax] = generate_beamforming_arrival_table(lats,deps,xducercoords2,nT,dY,dT,c0)

  idps_mat_ax=zeros(length(lats),length(deps));
  
  for ii=1:length(lats)
    lat=lats(ii);
    for jj=1:length(deps)
      dep=deps(jj);
      fcen=([lat/dY dep/dY]);
      [dd2 mdd2]=focusProfile2(fcen,xducercoords2,dT/dY*c0);
      idps_mat_ax(ii,jj)=mdd2;
    end
  end
  
