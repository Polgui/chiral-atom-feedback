clearvars

% vector: (ee;eg;ge:gg)

phiInput=(-pi:pi/100:pi);

Nt=100;
tspan=linspace(0,50,Nt+1);
dt=tspan(2)-tspan(1);

Delta=0;
DeltaL=1*Delta;
DeltaR=1*Delta;

gamma=1;
gammaL=1*gamma;
gammaR=1*gamma;

Omega=2;
OmegaL=1*Omega;

InitialState=[0,0,0,1];

rhoS0=InitialState'*InitialState;
rhoSvec0=mat2vec(rhoS0);

smL=[0,0,0,0;0,0,0,0;1,0,0,0;0,1,0,0];
smR=[0,0,0,0;1,0,0,0;0,0,0,0;0,0,1,0];
smT=sqrt(gammaL/(gammaL+gammaR))*smL + sqrt(gammaR/(gammaL+gammaR))*smR;

D=SOpre(smT)*SOpost(smT')-1/2*SOpre(smT'*smT)-1/2*SOpost(smT'*smT);
    
Purity_t=zeros(length(phiInput),Nt+1);
    
for iterphi=1:length(phiInput)
    
    OmegaR=1*Omega*exp(-1i*phiInput(iterphi));

    rhovect=zeros(length(rhoSvec0),Nt+1);
    rhovect(:,1)=rhoSvec0;

    Purity_t(iterphi,1)=trace(rhoS0*rhoS0);

    Hsys=-DeltaL*(smL'*smL)-DeltaR*(smR'*smR)-(OmegaL/2)*smL-((OmegaL/2)*smL)'...
        -(OmegaR/2)*smR-((OmegaR/2)*smR)' + 1i*sqrt(gammaL*gammaR/4)*(smL'*smR-smR'*smL);

    Liouvillian=-1i*SOpre(Hsys)+1i*SOpost(Hsys)+(gammaL+gammaR)*D;

    for k=1:Nt
        rhovect(:,k+1)=expv(dt,Liouvillian,rhovect(:,k));
        rho=vec2mat(rhovect(:,k+1));
        rho=(rho+rho')./2;
        rho=rho./trace(rho);
        Purity_t(iterphi,k+1)=trace(rho*rho);
    end

end
%%
figure('Position', [500, 200, 650, 400])
if length(phiInput)==1
    plot(tspan,Purity_t(:,1))
else
    [~, h2]=contourf(tspan,phiInput/pi,Purity_t,50);
    set(h2,'edgecolor','none')
    set(gca, 'fontsize',18)
    xlabel('t')
    ylabel('\phi')
    title(strcat('State purity, \gamma_L=',num2str(gammaL),', \gamma_R=',num2str(gammaR),...
        ', \Delta_L=',num2str(DeltaL),', \Delta_R=',num2str(DeltaR),...
        ', \Omega_L=',num2str(OmegaL),', \Omega_R=', num2str(abs(OmegaR)),'*exp(-i\phi)'))
    colormap('jet')
    colorbar
end