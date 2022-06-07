function [pxducer1 pxducer2]=fundharmpxducer(pxducer,f0,dT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: 2022-06-06
% LAST MODIFIED: 2022-06-06
% Filter fundamental and harmonic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nT=size(pxducer,1);

f=(0:nT-1)/(nT-1)/(dT);
%flt1=exp(-((f-f0)/f0*3).^4)';
flt1=exp(-((f-f0)/f0*4).^2)';
plot(f,flt1), xlim([0 4*f0]), grid on, hold on
%flt2=exp(-((f-2*f0)/f0*3).^4)';
flt2=exp(-((f-2*f0)/f0*4).^2)';

plot(f,flt2,'r'), xlim([0 4*f0]), grid on, hold off

flmat=flt1*ones(1,size(pxducer,2));
if(size(pxducer,3)==1)
  fpxducer=fft(double(pxducer));
  pxducer1=ifft(fpxducer.*flmat,'symmetric');
  pxducer1=single(pxducer1); 
end
if(size(pxducer,3)>1)
  for k=1:size(pxducer,3)
    fpxducer=fft(double(pxducer(:,:,k)));
    pxducer1(:,:,k)=single(ifft(fpxducer.*flmat,'symmetric'));
  end
end

flmat=flt2*ones(1,size(pxducer,2));
if(size(pxducer,3)==1)
  pxducer2=ifft(fpxducer.*flmat,'symmetric');
  pxducer2=single(pxducer2);
end
if(size(pxducer,3)>1)
  for k=1:size(pxducer,3)
    fpxducer=fft(double(pxducer(:,:,k)));
    pxducer2(:,:,k)=single(ifft(fpxducer.*flmat,'symmetric'));
  end
end
