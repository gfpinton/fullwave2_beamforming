%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% FIRST WRITTEN: 2018-06-21
% LAST MODIFIED: 2022-04-07
% Launch Fullwave 2 code, easy matlab wrapper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Basic variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c0=1540;         % average speed of sound (m/s)
rho0=1000;
f0=2e6;
omega0=2*pi*f0; % center radian frequency of transmitted wave
wY=8e-2;         % depth of simulation field (m)
duration=wY*2.3/c0;  % duration of simulation (s)
p0=5e5; % pressure in Pa
basedir='/kulm/scratch/fullwave2_beamforming_vishuman'
%%% Advanced variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfl=0.35;         % Courant-Friedrichs-Levi condition
ppw=12;  % number of points per spatial wavelength
dX=c0/omega0*2*pi/ppw % step size in x
dY=dX; % step size in y (please keep step sizes the same)
dT=dX/c0*cfl; % step size in time
xdc.nE=128; % number of transducer elements
xdc.ptch=6; % pitch in points
xdc.kerf=0; %kerf in points
xdc.width=xdc.ptch*dX*xdc.nE; % (m)
ppp=ppw/cfl;     % points per period
wX=xdc.width+1e-3;         % width of simulation field (m)
%%% Grid size calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda=c0/omega0*2*pi; % wavelength (m)
nX=round(wX/lambda*ppw);  % number of lateral elements
nY=round(wY/lambda*ppw);  % number of depth elements
nT=round(duration*c0/lambda*ppw/cfl); % number of time points
modT=3; modX=4; modY=4;
dT2=dT*modT;
nT2=length(1:modT:nT); nX2=length(1:modX:nX); nY2=length(1:modY:nY);
%%MAPS AND GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mpg.cfl=cfl;
mpg.dX=dX; mpg.dY=dY; mpg.dT=dT;
mpg.nX=nX; mpg.nY=nY; mpg.nT=nT;
mpg.modX=modX;mpg.modY=modY;mpg.modT=modT;
mpg.nX2=nX2; mpg.nY2=nY2; mpg.nT=nT2;
mpg.c=ones(nX,nY)*1540;   % speed of sound map (m/s)
mpg.rho=ones(nX,nY)*1000; % density map (kg/m^3)
mpg.A=ones(nX,nY)*0.5;    % attenuation map (dB/MHz/cm)
mpg.bovera=7.6*ones(nX,nY);% -2*ones(nX,nY);    % nonlinearity map
load vishuman_abdominal_slice
figure(1), imagesc(cmap'), figure(2), imagesc(rhomap')
cmap=interp2easy(cmap,dm/dX,dm/dX,'nearest');
rhomap=interp2easy(rhomap,dm/dX,dm/dX,'nearest');
Amap=interp2easy(Amap,dm/dX,dm/dX,'nearest');
betamap=interp2easy(Nmap,dm/dX,dm/dX,'nearest');

cmap=cmap(round(end/2)-round(nX/2)+1:round(end/2)-round(nX/2)+nX,1:nY);
rhomap=rhomap(round(end/2)-round(nX/2)+1:round(end/2)-round(nX/2)+nX,1:nY);
Amap=Amap(round(end/2)-round(nX/2)+1:round(end/2)-round(nX/2)+nX,1:nY);
betamap=betamap(round(end/2)-round(nX/2)+1:round(end/2)-round(nX/2)+nX,1:nY);
gfilt=(5/10)^2*ppw/2; % correct for pixelization
cmap=imgaussfilt(cmap,gfilt);
rhomap=imgaussfilt(rhomap,gfilt);
Amap=imgaussfilt(Amap,gfilt);
betamap=imgaussfilt(betamap,gfilt);
mpg.c=cmap; mpg.A=Amap; mpg.betamap=betamap;

% scatterers
mpg.scat_density=0.15;
mpg.scats=rand(nX,nY);
mpg.scats(find(mpg.scats>mpg.scat_density))=0; mpg.scats=mpg.scats/max(max(mpg.scats));
mean(mean(mpg.scats))
mpg.scats(:,1:10)=0; % don't put scatters inside your transducer
mpg.rhosr=0.0375*2; % scatterer impedance contrast
for k=1:10
  idl=circleIdx(size(mpg.scats),[nX/2 k*1e-2/dX],3e-3/dY);
  mpg.scats(idl)=0; % anechoic region
end
for k=1:10
  idl=circleIdx(size(mpg.scats),[nX/2-5e-3/dX k*1e-2/dX],0.65);
  mpg.scats(idl)=0.8/mpg.rhosr; % anechoic region
end
for k=1:0.5:10
  idl=circleIdx(size(mpg.scats),[nX/2+5e-3/dX k*1e-2/dX],0.65);
  mpg.scats(idl)=0.8/mpg.rhosr; % anechoic region
end
imagesc(mpg.scats')

mpg.rho=rhomap-mpg.scats*1000*mpg.rhosr;
imagesc(mpg.rho'), axis equal, colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define transmit characteristics %%
tx.nTx=40; % number of Tx events
tx.dep=6e-2/dX; % focal depth (pixels)
tx.fnumber=tx.dep/xdc.width*dX;
tx.bmw=floor(tx.fnumber*lambda/dX/2)% lambda*Z/2D beamwidth (pixels)
tx.bmw=lambda/2/dX; % beamwidth (pixels), is this an integer multiple of dX?
if(mod(tx.bmw,1)), disp('Warning, tx.bmw not an integer multiple of dX'), end
tx.apex=[nX/2 3];
tx.ncycles=2; % number of cycles in pulse
tx.dur=1; % exponential drop-off of envelope

%%% Generate input coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%
inmap=zeros(nX,nY); 
inmap(:,1:8)=ones(nX,8);
imagesc(inmap'), axis equal, axis tight
incoords=mapToCoords(inmap); % note zero indexing for compiled code
plot(incoords(:,1),incoords(:,2),'.')
%% COORDINATES INPUT OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xdc.incoords xdc.outcoords xdc.incoords2 xdc.outcoords2]=genxdccoords(xdc.nE,xdc.ptch,nX,xdc.kerf);
plot(xdc.incoords(:,1),xdc.incoords(:,2),'.'), hold on
plot(xdc.incoords2(:,1),xdc.incoords2(:,2),'.')
plot(xdc.outcoords(:,1),xdc.outcoords(:,2),'.')
plot(xdc.outcoords2(:,1),xdc.outcoords2(:,2),'.'), hold off
xdc.outcoords_field = coordsMatrix(mpg.nX,mpg.nY,mpg.modX,mpg.modY); xdc.outcoords_field(:,3)=-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[tx.focs ]=genlinearfocs(tx.nTx,tx.bmw,tx.dep,tx.apex)
plot(xdc.incoords(:,1),xdc.incoords(:,2),'.'),hold on
plot(xdc.outcoords(:,1),xdc.outcoords(:,2),'.')
plot(tx.focs(:,1),tx.focs(:,2),'.'),
plot(tx.apex(1),tx.apex(2),'.')
hold off, axis equal, grid on
legend('incoords','outcoords','focs','apex')
%% generate focal delays %%
tx.dd=focusDelays(tx.focs,xdc.incoords,cfl,0); tx.dd=tx.dd-min(min(tx.dd));
for tt=1:xdc.nE % one focal delay per element
  idt=find(xdc.incoords(:,3)==tt);
  tx.dd(idt,:)=ones(length(idt),1)*mean(tx.dd(idt,:));
end
tx.nTic=round(max(max(tx.dd))+tx.ncycles*5*ppw/cfl)
t=(0:nT-1)/nT*duration-tx.ncycles/omega0*2*pi;
icvec=exp(-(1.05*t*omega0/(tx.ncycles*pi)).^(2*tx.dur)).*sin(t*omega0)*p0;
tx.icvec=icvec(1:tx.nTic);

%% REFERENCE SIM %%
outdir=[ basedir '/']
eval(['!mkdir -p ' outdir]);
n=1;
icmat=focusCoordsDD(round(tx.dd(:,n)),tx.icvec,0);
icmat(find(xdc.incoords(:,3)==0),:)=0; % zero out elements labeled zero
%%% need to implement transmit fnumber here %%%
%%% element averaging on transmit %%

cwd=pwd; addpath(cwd);
cd(outdir)
launch_fullwave2_try6_nln_relaxing4(c0,omega0,wX,wY,duration,p0,ppw,cfl,mpg.c,mpg.rho,mpg.A,1+mpg.betamap,xdc.incoords,xdc.outcoords,icmat); writeVabs('int',modT,'modT');
eval(['!rm icmat.dat'])
cd(cwd);

%% INDIVIDUAL TRANSMIT-RECEIVE SIMS %%%%%%%%%%%%%%
for n=1:tx.nTx
  outdir=[ basedir '/txrx_' num2str(n)]
  eval(['!mkdir -p ' outdir]);
  eval(['!cp fullwave2_try6_nln_relaxing ' outdir]);
  
  icmat=focusCoordsDD(round(tx.dd(:,n)),tx.icvec,0);
  icmat(find(xdc.incoords(:,3)==0),:)=0; % zero out elements labeled zero

  %for tt=1:xdc.nE
  %  idc=find(xdc.outcoords(:,3)==tt);
  %end
 
  
   %idx=find(abs(xdc.incoords_thetas-tx.thetas(n))>tx.fthetanumber/2);
   %icmat(idx,:)=0;
%idx=find((abs(xdc.incoords(:,1)-tx.focs(n,1)))/tx.focs(n,2)>tx.fnumber/2);
%icmat(idx,:)=0;
  imagesc(icmat), drawnow
  
  cwd=pwd; addpath(cwd);
  cd(outdir)
  eval(['!ln -s ../*.dat .'])
  writeIC('icmat.dat',(icmat)');
  %eval('!./fullwave2_try6_nln_relaxing & ') %% !! UNCOMMENT TO RUN EXECUTABLE
  cd(cwd);
end

%% read pxducer %%%
idc=find(xdc.outcoords(:,3)>=1); length(idc)
pxducer=zeros(nT2,length(idc),tx.nTx,'single');
pxducer_avg=zeros(nT2,xdc.nE,tx.nTx,'single');
 for ii=1:tx.nTx
   outdir=[basedir '/txrx_' num2str(ii) '/']
   nRun=sizeOfFile([outdir 'genout.dat'])/4/size(xdc.outcoords,1)
   if(nRun>0)
     tmp=readGenoutSlice([outdir 'genout.dat'],0:nRun-1,size(xdc.outcoords,1),idc);
     tmp(end:nT2,:)=0;
     if(nRun>nT2)
       tmp=tmp(1:nT2,:)*0; disp('WARNING ZEROING PXDUCER DUE TO OVERFLOW ERROR')
     end
     pxducer(:,:,ii)=tmp;     
   end
 end
 
 for ii=1:tx.nTx
   for tt=1:xdc.nE
     tmp=pxducer(:,find(xdc.outcoords(:,3)==tt),ii);
     tmp=mean(tmp,2);
%pxducer_avg(:,tt,ii)=mean(pxducer(:,find(xdc.outcoords(:,4)==tt)),2);
     pxducer_avg(:,tt,ii)=tmp;
   end
 end
 
 imagesc(powcompress(pxducer(:,:,round(end/2)),1/3))
 imagesc(powcompress(pxducer_avg(:,:,round(end/2)),1/3))

 %[pxducer1 pxducer2]=fundharmpxducer(pxducer_avg,f0,dT*modT);
 
%%% GENERATE IMAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 deps=1e-3:lambda/8:nY*dY/1.1;
 lats = tx.focs(:,1)*dX;
 [val idt0s] = max(abs(hilbert(pxducer)));
 idt0s=squeeze(idt0s); imagesc(idt0s)
 idt0s=squeeze(max(idt0s))
 
 
 idx=find(xdc.outcoords2(:,1));
 rx.idps_mat_rx=zeros(length(lats),length(deps),size(xdc.outcoords2,1)+1);
 for ii=1:length(lats)
   for jj=1:length(deps)
     fcen=[lats(ii) deps(jj)]/dX;
     [dd mdd]=focusProfile2(fcen,xdc.outcoords2(idx,:),dT/dY*c0*modT);
     idp=double((nT2*(idx-1))+mdd); %no rounding
     rx.idps_mat_rx(ii,jj,1)=length(idp);
     rx.idps_mat_rx(ii,jj,2:length(idp)+1)=(idp);
   end
 end
 
 rx.idps_mat_ax=zeros(length(lats),length(deps));
 for ii=1:length(lats)
   for jj=1:length(deps)
     fcen=[lats(ii) deps(jj)]/dX;
     [dd2 mdd2]=focusProfile2(fcen,tx.apex,dT/dY*c0*modT);
     rx.idps_mat_ax(ii,jj)=mdd2;
   end
 end
 rx.idps_mat_ax=rx.idps_mat_ax-min(min(rx.idps_mat_ax));
 
 
 %% delay and sum %%
 bm=zeros(length(lats),length(deps));
 bm2=zeros(length(lats),length(deps));
 for ii=1:length(lats)
   pxducer_now=pxducer_avg(:,:,ii);
       %  px=pxducer(:,round(end/2)); [val idt0]=max(abs(hilbert(px)))
   for jj=1:length(deps)
     bm(ii,jj)=sum(pxducer_now(round(rx.idps_mat_rx(ii,jj,2:rx.idps_mat_rx(ii,jj,1))+rx.idps_mat_ax(ii,jj)+idt0s(ii))));
   end
   imagesc(dbah(bm(:,:)'),[-60 0]);
   colormap gray, title(num2str(n)), drawnow
end

%% PLOT THE BMODE IMAGE %%
figure(1)
subplot(1,2,1)
imagesc(lats*1e3,deps*1e3,dbah(squeeze(bm')),[-45 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
subplot(1,2,2)
imagesc(lats*1e3,deps*1e3,dbah(squeeze(bm')),[-45 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
axis([mean(lats*1e3) mean(lats*1e3)+6 55 65])
title('Reference')

interpfac=8;
figure(2)
subplot(1,2,1)
imagesc(lats*1e3,deps*1e3,dbah(interp2easy(bm',interpfac,interpfac)),[-45 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
subplot(1,2,2)
imagesc(lats*1e3,deps*1e3,dbah(interp2easy(bm',interpfac,interpfac)),[-45 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
axis([mean(lats*1e3) mean(lats*1e3)+6 55 65])
title('Interpolation of B-mode')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENVELOPE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
subplot(1,2,1)
imagesc(lats*1e3,deps*1e3,dbenv(interp2easy(bm',interpfac,interpfac)),[-45 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
subplot(1,2,2)
imagesc(lats*1e3,deps*1e3,dbenv(interp2easy(bm',interpfac,interpfac)),[-45 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
axis([mean(lats*1e3) mean(lats*1e3)+6 55 65])
title('60-tap Hilbert envelope')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPRESSION SCALE AND METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
subplot(1,2,1)
imagesc(lats*1e3,deps*1e3,dbah(interp2easy(bm',interpfac,interpfac)),[-35 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
subplot(1,2,2)
imagesc(lats*1e3,deps*1e3,dbah(interp2easy(bm',interpfac,interpfac)),[-35 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
axis([mean(lats*1e3) mean(lats*1e3)+6 55 65])
title('alternate dB Scale')

figure(5)
subplot(1,2,1)
imagesc(lats*1e3,deps*1e3,dbah(interp2easy(bm',interpfac,interpfac)),[-45 -10])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
subplot(1,2,2)
imagesc(lats*1e3,deps*1e3,dbah(interp2easy(bm',interpfac,interpfac)),[-45 -10])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
axis([mean(lats*1e3) mean(lats*1e3)+6 55 65])
title('alternate dB Scale')

figure(6)
subplot(1,2,1)
imagesc(lats*1e3,deps*1e3,pwah(interp2easy(bm',interpfac,interpfac),1/3),[-70 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
subplot(1,2,2)
imagesc(lats*1e3,deps*1e3,pwah(interp2easy(bm',interpfac,interpfac),1/3),[-70 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
axis([mean(lats*1e3) mean(lats*1e3)+6 55 65])
title('Power 1/3 Scale')

figure(7)
subplot(1,2,1)
imagesc(lats*1e3,deps*1e3,pwah(interp2easy(bm',interpfac,interpfac),1/4),[-20 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
subplot(1,2,2)
imagesc(lats*1e3,deps*1e3,pwah(interp2easy(bm',interpfac,interpfac),1/4),[-20 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
axis([mean(lats*1e3) mean(lats*1e3)+6 55 65])
title('Power 1/4 Scale')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENVELOPE BEFORE DETECTING      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% delay and sum %%
bm_ebd=zeros(length(lats),length(deps));
for ii=1:length(lats)
  pxducer_now=pxducer_avg(:,:,ii);
  pxducer_now=hilbert(pxducer_now);
       %  px=pxducer(:,round(end/2)); [val idt0]=max(abs(hilbert(px)))
  for jj=1:length(deps)
     bm_ebd(ii,jj)=sum(pxducer_now(round(rx.idps_mat_rx(ii,jj,2:rx.idps_mat_rx(ii,jj,1))+rx.idps_mat_ax(ii,jj)+idt0s(ii))));
  end
  %imagesc(bm(:,:)',[-45 0]);
  %colormap gray, title(num2str(n)), drawnow
end
bm_ebd=abs(bm_ebd); 

interpfac=8;
figure(8)
subplot(1,2,1)
imagesc(lats*1e3,deps*1e3,interp2easy(dbzero(bm_ebd'),interpfac,interpfac),[-45 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
subplot(1,2,2)
imagesc(lats*1e3,deps*1e3,interp2easy(dbzero(bm_ebd'),interpfac,interpfac),[-45 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
axis([mean(lats*1e3) mean(lats*1e3)+6 55 65])
title('Hilbert before detecting')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEAMFORMING SPEED OF SOUND    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cbm=mean(mean(mpg.c))

idx=find(xdc.outcoords2(:,1));
rx.idps_mat_rx=zeros(length(lats),length(deps),size(xdc.outcoords2,1)+1);
for ii=1:length(lats)
  for jj=1:length(deps)
    fcen=[lats(ii) deps(jj)]/dX;
    [dd mdd]=focusProfile2(fcen,xdc.outcoords2(idx,:),dT/dY*cbm*modT);
    idp=double((nT2*(idx-1))+mdd); %no rounding
    rx.idps_mat_rx(ii,jj,1)=length(idp);
    rx.idps_mat_rx(ii,jj,2:length(idp)+1)=(idp);
  end
end

rx.idps_mat_ax=zeros(length(lats),length(deps));
for ii=1:length(lats)
  for jj=1:length(deps)
    fcen=[lats(ii) deps(jj)]/dX;
    [dd2 mdd2]=focusProfile2(fcen,tx.apex,dT/dY*cbm*modT);
    rx.idps_mat_ax(ii,jj)=mdd2;
  end
end
rx.idps_mat_ax=rx.idps_mat_ax-min(min(rx.idps_mat_ax));


 %% delay and sum %%
bm=zeros(length(lats),length(deps));
for ii=1:length(lats)
  pxducer_now=pxducer_avg(:,:,ii);
       %  px=pxducer(:,round(end/2)); [val idt0]=max(abs(hilbert(px)))
  for jj=1:length(deps)
     bm(ii,jj)=sum(pxducer_now(round(rx.idps_mat_rx(ii,jj,2:rx.idps_mat_rx(ii,jj,1))+rx.idps_mat_ax(ii,jj)+idt0s(ii))));
  end
end

figure(9)
subplot(1,2,1)
imagesc(lats*1e3,deps*1e3,dbah(interp2easy(bm',interpfac,interpfac)),[-45 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
subplot(1,2,2)
imagesc(lats*1e3,deps*1e3,dbah(interp2easy(bm',interpfac,interpfac)),[-45 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
axis([mean(lats*1e3) mean(lats*1e3)+6 55 65])
title('Mean speed of sound ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cbm=c0-c0*0.05

idx=find(xdc.outcoords2(:,1));
rx.idps_mat_rx=zeros(length(lats),length(deps),size(xdc.outcoords2,1)+1);
for ii=1:length(lats)
  for jj=1:length(deps)
    fcen=[lats(ii) deps(jj)]/dX;
    [dd mdd]=focusProfile2(fcen,xdc.outcoords2(idx,:),dT/dY*cbm*modT);
    idp=double((nT2*(idx-1))+mdd); %no rounding
    rx.idps_mat_rx(ii,jj,1)=length(idp);
    rx.idps_mat_rx(ii,jj,2:length(idp)+1)=(idp);
  end
end

rx.idps_mat_ax=zeros(length(lats),length(deps));
for ii=1:length(lats)
  for jj=1:length(deps)
    fcen=[lats(ii) deps(jj)]/dX;
    [dd2 mdd2]=focusProfile2(fcen,tx.apex,dT/dY*cbm*modT);
    rx.idps_mat_ax(ii,jj)=mdd2;
  end
end
rx.idps_mat_ax=rx.idps_mat_ax-min(min(rx.idps_mat_ax));


 %% delay and sum %%
bm=zeros(length(lats),length(deps));
for ii=1:length(lats)
  pxducer_now=pxducer_avg(:,:,ii);
       %  px=pxducer(:,round(end/2)); [val idt0]=max(abs(hilbert(px)))
  for jj=1:length(deps)
     bm(ii,jj)=sum(pxducer_now(round(rx.idps_mat_rx(ii,jj,2:rx.idps_mat_rx(ii,jj,1))+rx.idps_mat_ax(ii,jj)+idt0s(ii))));
  end
end

figure(10) 
subplot(1,2,1)
imagesc(lats*1e3,deps*1e3,dbah(interp2easy(bm',interpfac,interpfac)),[-45 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
subplot(1,2,2)
imagesc(lats*1e3,deps*1e3,dbah(interp2easy(bm',interpfac,interpfac)),[-45 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
axis([mean(lats*1e3) mean(lats*1e3)+6 55 65])
title('Mean speed of sound -5% ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AXIAL UNDERSAMPLING   %%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modT2=5*modT;
pxducer_avg2=pxducer_avg(1:modT2/modT:end,:,:);

idx=find(xdc.outcoords2(:,1));
rx.idps_mat_rx=zeros(length(lats),length(deps),size(xdc.outcoords2,1)+1);
for ii=1:length(lats)
  for jj=1:length(deps)
    fcen=[lats(ii) deps(jj)]/dX;
    [dd mdd]=focusProfile2(fcen,xdc.outcoords2(idx,:),dT/dY*c0*modT2);
    idp=double((size(pxducer_avg2,1)*(idx-1))+mdd); %no rounding
    rx.idps_mat_rx(ii,jj,1)=length(idp);
    rx.idps_mat_rx(ii,jj,2:length(idp)+1)=(idp);
  end
end

rx.idps_mat_ax=zeros(length(lats),length(deps));
for ii=1:length(lats)
  for jj=1:length(deps)
    fcen=[lats(ii) deps(jj)]/dX;
    [dd2 mdd2]=focusProfile2(fcen,tx.apex,dT/dY*c0*modT2);
    rx.idps_mat_ax(ii,jj)=mdd2;
  end
end
rx.idps_mat_ax=rx.idps_mat_ax-min(min(rx.idps_mat_ax));


 %% delay and sum %%
bm=zeros(length(lats),length(deps));
for ii=1:length(lats)
  pxducer_now=pxducer_avg2(:,:,ii);
       %  px=pxducer(:,round(end/2)); [val idt0]=max(abs(hilbert(px)))
  for jj=1:length(deps)
     bm(ii,jj)=sum(pxducer_now(round(rx.idps_mat_rx(ii,jj,2:rx.idps_mat_rx(ii,jj,1))+rx.idps_mat_ax(ii,jj)+idt0s(ii)*modT/modT2)));
  end
end

figure(11)
subplot(1,2,1)
imagesc(lats*1e3,deps*1e3,dbah(interp2easy(bm',interpfac,interpfac)),[-45 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
subplot(1,2,2)
imagesc(lats*1e3,deps*1e3,dbah(interp2easy(bm',interpfac,interpfac)),[-45 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
axis([mean(lats*1e3) mean(lats*1e3)+6 55 65])
title('Axial undersampling')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% idt0 ERRORS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idt0s_err=rand(size(idt0s))*ppp/modT*2;idt0s_err=idt0s_err-mean(idt0s_err);
plot(idt0s), hold on
plot(idt0s+idt0s_err), hold off, grid on
legend('idt0','idt0+error'), xlabel('Transmit event'), ylabel('idt0 (time pixel)')

idx=find(xdc.outcoords2(:,1));
rx.idps_mat_rx=zeros(length(lats),length(deps),size(xdc.outcoords2,1)+1);
for ii=1:length(lats)
  for jj=1:length(deps)
    fcen=[lats(ii) deps(jj)]/dX;
    [dd mdd]=focusProfile2(fcen,xdc.outcoords2(idx,:),dT/dY*c0*modT);
    idp=double((nT2*(idx-1))+mdd); %no rounding
    rx.idps_mat_rx(ii,jj,1)=length(idp);
    rx.idps_mat_rx(ii,jj,2:length(idp)+1)=(idp);
  end
end

rx.idps_mat_ax=zeros(length(lats),length(deps));
for ii=1:length(lats)
  for jj=1:length(deps)
    fcen=[lats(ii) deps(jj)]/dX;
    [dd2 mdd2]=focusProfile2(fcen,tx.apex,dT/dY*c0*modT);
    rx.idps_mat_ax(ii,jj)=mdd2;
  end
end
rx.idps_mat_ax=rx.idps_mat_ax-min(min(rx.idps_mat_ax));


 %% delay and sum %%
bm=zeros(length(lats),length(deps));
for ii=1:length(lats)
  pxducer_now=pxducer_avg(:,:,ii);
       %  px=pxducer(:,round(end/2)); [val idt0]=max(abs(hilbert(px)))
  for jj=1:length(deps)
     bm(ii,jj)=sum(pxducer_now(round(rx.idps_mat_rx(ii,jj,2:rx.idps_mat_rx(ii,jj,1))+rx.idps_mat_ax(ii,jj)+idt0s(ii)+idt0s_err(ii))));
  end
end

figure(12)
subplot(1,2,1)
imagesc(lats*1e3,deps*1e3,dbah(interp2easy(bm',interpfac,interpfac)),[-45 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
subplot(1,2,2)
imagesc(lats*1e3,deps*1e3,dbah(interp2easy(bm',interpfac,interpfac)),[-45 0])
colormap gray, colorbar
xlabel('mm'), ylabel('mm')
axis equal, axis tight
axis([mean(lats*1e3) mean(lats*1e3)+6 55 65])
title('idt0 error')

