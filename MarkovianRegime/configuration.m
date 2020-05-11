clearvars
addpath('functions')
% vector: (l;r;g)

phiInput=(pi*0.98:pi/2000:pi);
OmegaInput=(0.1:0.1:60);

N_trial=1;

Delta=0;
DeltaL=1*Delta;
DeltaR=1*Delta;

gamma=1;
gammaL=1.*gamma;
gammaR=1.*gamma;

rangeimpL=0.;
rangeimpR=rangeimpL;

gammap=gamma*(1/0.9999 - 1);

phiRightCoupling=0.*pi;
   

Purity_t=zeros(length(phiInput),length(OmegaInput));
excitationL=zeros(length(phiInput),length(OmegaInput));
excitationR=zeros(length(phiInput),length(OmegaInput));
excitationG=zeros(length(phiInput),length(OmegaInput));
    
for iterphi=1:length(phiInput)
   for iterOmega=1:length(OmegaInput)
       
       phiPhotons=phiInput(iterphi);
       
        OmegaL=OmegaInput(iterOmega);
        OmegaR=OmegaL;
        rho=zeros(3,3);
        for itertrial=1:N_trial
            
        smL=[0,0,0;0,0,0;1,0,0];
        smR=[0,0,0;0,0,0;0,1,0];

        impL=rangeimpL;
        impR=rangeimpR;
        
        smbL=sqrt(gammaL*(1-impL))*smL+sqrt(gammaR*impR)*smR*exp(1i*phiRightCoupling);
        smbR=sqrt(gammaR*(1-impR))*smR*exp(1i*phiRightCoupling)+sqrt(gammaL*impL)*smL;

        DL=SOpre(smL)*SOpost(smL')-1/2*SOpre(smL'*smL)-1/2*SOpost(smL'*smL);
        DR=SOpre(smR)*SOpost(smR')-1/2*SOpre(smR'*smR)-1/2*SOpost(smR'*smR);
        
        smbT=smbL+exp(1i*phiPhotons)*smbR;
        DbT=SOpre(smbT)*SOpost(smbT')-1/2*SOpre(smbT'*smbT)-1/2*SOpost(smbT'*smbT);

        Hsys=-DeltaL*(smL'*smL)-DeltaR*(smR'*smR)-(OmegaL/2)*smL-((OmegaL/2)*smL)'...
            -(OmegaR/2)*smR-((OmegaR/2)*smR)'; 
        
        Hdd=(1i/2)*(exp(1i*phiPhotons)*smbL'*smbR - exp(-1i*phiPhotons)*smbR'*smbL);

        Liouvillian=-1i*SOpre(Hsys+Hdd)+1i*SOpost(Hsys+Hdd)+gammap*DL+gammap*DR+DbT;

        [U,V]=eig(Liouvillian,'vector');
        if length(V(abs(V)<1e-6))~= 1
            error('The Liouvillian has a number of 0-eigenvalues different than 1')
        end
        if ~isempty(V(V>1e-4))
            error('The Liouvillian has positive eigenvalues')
        end   
        rhovect=U(:,abs(V)<1e-6);

% rhovect(:)=expm(Liouvillian*1000)*[0 0 0 0 0 0 0 0 1]';

        rho_trial=vec2mat(rhovect(:));
        rho=rho+(rho_trial./trace(rho_trial))/N_trial;
        end
        
        Purity_t(iterphi,iterOmega)=trace(rho*rho);
        excitationL(iterphi,iterOmega)=rho(1,1);
        excitationR(iterphi,iterOmega)=rho(2,2);
        excitationG(iterphi,iterOmega)=rho(3,3);
   end
end
%%
figure('Position', [500, 200, 650, 400])

[~, h2]=contourf(OmegaInput,phiInput/pi,real(Purity_t),800);
caxis([0.6 1])
set(h2,'edgecolor','none')
set(gca, 'fontsize',18)
xlabel('\Omega')
ylabel('\phi')
title(strcat('State purity, \gamma_L=',num2str(gammaL),', \gamma_R=',num2str(gammaR),...
    ', \Delta_L=',num2str(DeltaL),', \Delta_R=',num2str(DeltaR),...
    ', \Omega_L=\Omega, \Omega_R=\Omega*exp(i\phi)'))
colormap('jet')
colorbar

% figure('Position', [600, 200, 650, 400])
% 
% [~, h2]=contourf(OmegaInput,phiInput/pi,real(excitationL),800);
% set(h2,'edgecolor','none')
% set(gca, 'fontsize',18)
% xlabel('\Omega')
% ylabel('\phi')
% title(strcat('Excitation L, \gamma_L=',num2str(gammaL),', \gamma_R=',num2str(gammaR),...
%     ', \Delta_L=',num2str(DeltaL),', \Delta_R=',num2str(DeltaR),...
%     ', \Omega_L=\Omega, \Omega_R=\Omega*exp(i\phi)'))
% colormap('jet')
% colorbar
% 
% figure('Position', [700, 200, 650, 400])
% 
% [~, h2]=contourf(OmegaInput,phiInput/pi,real(excitationR),800);
% set(h2,'edgecolor','none')
% set(gca, 'fontsize',18)
% xlabel('\Omega')
% ylabel('\phi')
% title(strcat('Excitation R, \gamma_L=',num2str(gammaL),', \gamma_R=',num2str(gammaR),...
%     ', \Delta_L=',num2str(DeltaL),', \Delta_R=',num2str(DeltaR),...
%     ', \Omega_L=\Omega, \Omega_R=\Omega*exp(i\phi)'))
% colormap('jet')
% colorbar
% 
% figure('Position', [800, 200, 650, 400])
% 
% [~, h2]=contourf(OmegaInput,phiInput/pi,real(excitationG),800);
% set(h2,'edgecolor','none')
% set(gca, 'fontsize',18)
% xlabel('\Omega')
% ylabel('\phi')
% title(strcat('Excitation G, \gamma_L=',num2str(gammaL),', \gamma_R=',num2str(gammaR),...
%     ', \Delta_L=',num2str(DeltaL),', \Delta_R=',num2str(DeltaR),...
%     ', \Omega_L=\Omega, \Omega_R=\Omega*exp(i\phi)'))
% colormap('jet')
% colorbar