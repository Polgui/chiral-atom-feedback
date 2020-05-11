%% SETTINGS
clearvars
addpath('QSSE_functions/analysis/')
addpath('QSSE_functions/evolution/')
addpath('QSSE_functions/others/')
addpath('QSSE_functions/unitaries/')
 
scriptname='NoCavity_Omega=2';
disp(scriptname)
disp('CURRENT PID:')
feature getpid
 
c2 = parcluster('local');
c2.NumWorkers = 12;
mypool=parpool(c2.NumWorkers); 
  
Totaltime=100;        % Total simulation time
dt=0.02;           % Timestep

TauMin = dt;
dtau = dt;
TauMax = 0.1;
Tau = (TauMin:dtau:TauMax);

CorrelationTime=2*TauMax;
    
phiMin=-pi;
dphi=pi/8;
phiMax=pi;
phiInput=(phiMin:dphi:phiMax);

N_wav=1;                                % Number of waveguides
N_sys=1;                                % Number of separated components in the system
e_init=[1 0 0];                         % Initial excitation of each state

dimB=3;                                 % Dimension for bosonic modes = 3 (no photon, 1 photon, 2 photons) for each waveguide
dimS=3;                                 % Dimension for each atom

Nt=round(Totaltime/dt);    % Total number of timesteps

maxSchmidtrank=256; % Maximum Schmidt rank for the SVD 
Eps = 10^(-5);      % Singular values truncation   

% Decay rates in the waveguide
gamma=1.;
gammaL=1*gamma; 
gammaR=1*gamma;

OmegaL=2;

DeltaL=1;
DeltaR=1;
    
showprogress=0;
showfigures=1;
saveresults=1;

excitations=cell(1,length(phiInput));
Lambdaloop=cell(1,length(phiInput));
S_profile=cell(1,length(phiInput));
auto_corr_1=cell(1,length(phiInput));
auto_corr_2=cell(1,length(phiInput));
Purity=cell(1,length(phiInput));
PhotonsOut=cell(1,length(phiInput));
PhotonsIn=cell(1,length(phiInput));
PhotonNumberLoop=cell(1,length(phiInput));

for iterphi=1:length(phiInput)
display(['iterphi:' num2str(iterphi)])    
par_excitations=cell(1,length(Tau));
par_Lambdaloop=cell(1,length(Tau));
par_S_profile=cell(1,length(Tau));
par_auto_corr_1=cell(1,length(Tau));
par_auto_corr_2=cell(1,length(Tau));
par_Purity=cell(1,length(Tau));
par_PhotonsOut=cell(1,length(Tau));
par_PhotonsIn=cell(1,length(Tau));
par_PhotonNumberLoop=cell(1,length(Tau));

myphi=phiInput(iterphi);

parfor iterTau=1:length(Tau)
    tau=Tau(iterTau);              % Delay time tau=d/v
    display(['iterTau:' num2str(iterTau)])
    m=cell(1,1);
    Theta=cell(1,1);
    UMPS=cell(1,1);
    
    % Number of steps between the interacting components
    m{1}=[-max(ceil(tau/dt),1) +1];

    % Losses associated
    Theta{1}=[0 0];

    %% UNITARY

    OmegaR=OmegaL*exp(1i*myphi);

    % Definition of Unitary

    UMPS{1}=V_system([gammaL gammaR], [OmegaL OmegaR], [DeltaL DeltaR], dimB,dt);

    %% INITIAL STATE

    m_max=1;
    N_loss=0;
    for j=1:length(m)
        m_max=max(m_max,max(max(abs(m{j}))));
        N_loss=N_loss+length(Theta{j}(Theta{j}~=0));
    end
    m_max=(N_wav+N_loss)*m_max;

    psiB=zeros(1,dimB);psiB(1)=1;           % Initial state of bosonic modes: vacuum
    Gammavac=reshape(psiB,[1,1,dimB]);
    Gammain=cell(1,N_wav);
    for j=1:N_wav+N_loss
        Gammain{j}=Gammavac;
    end

    GammaS=cell(1,N_sys);
    % Initial system state. [0,1]:excited
    for j=1:N_sys
        GammaS{j}=reshape(e_init,[1,1,dimS]);
    end


    Gamma0=cell(1,m_max);   % Initial Gamma
    Lambda0=cell(1,m_max);  % Initial Lambda
    for n=1:m_max+1
        Gamma0{n}=Gammavac;
        Lambda0{n}=1;
    end

    for j=1:N_sys
        Gamma0{m_max+j}=GammaS{j};
        Lambda0{m_max+j}=1;
    end

    %% SIMULATION

    [Gamma,Lambda,rhoS,rhoField, Error_Tracker, LambdaLoop]...
        = FullEvolution(Gamma0,Lambda0,Gammain,UMPS,m, N_wav,Theta,N_loss, Nt,maxSchmidtrank, Eps, showprogress);

    %% ANALYSIS

    % Evolution of the system populations
    exc_ev=zeros(3,Nt+1); 
    exc_ev(:,1)=abs(e_init).^2;
    [exc_ev]=GetExcitationProbabilities(rhoS,exc_ev);
    par_excitations{iterTau}=exc_ev;
    
    % Evolution of the system purity
    SystemPurity = zeros(1,Nt);
    for i=1:Nt
        SystemPurity(i) = trace(rhoS{1,i}*rhoS{1,i}); 
    end
    par_Purity{iterTau}=SystemPurity;
    
    % Evolution of the mean number of photons emitted during dt
    PhotonEmission = zeros(1,Nt);
    PhotonN = zeros(1,dimB-1);
    for i=1:Nt
        for j=1:dimB-1
            PhotonN(j) = rhoField{1,1,i}(j+1,j+1);
        end
        PhotonEmission(i) = PhotonN * (1:dimB-1)';
    end
    par_PhotonsOut{iterTau} = PhotonEmission;
    
    
    % Evolution of the mean number of photons injected in the loop during dt
    PhotonInjection = zeros(1,Nt);
    PhotonN = zeros(1,dimB-1);
    
    for i=1:Nt
        for j=1:dimB-1
            PhotonN(j) = rhoField{1,2,i}(j+1,j+1);
        end
        PhotonInjection(i) = PhotonN * (1:dimB-1)';
    end
    
    par_PhotonsIn{iterTau}=PhotonInjection;
    

    % Evolution of the number of photons in the loop
    N_loop=zeros(1,Nt+1);
    for i=1:Nt
        N_loop(i+1)=N_loop(i)+PhotonInjection(i);
        if (i-ceil(tau/dt)>0)
           N_loop(i+1)=N_loop(i+1)-PhotonInjection(i-ceil(tau/dt)); 
        end
    end
    par_PhotonNumberLoop{iterTau}=N_loop;
    
    % Evolution of the entropy profile in the loop
    par_Lambdaloop{iterTau} = LambdaLoop;

    % Entropy profile of single cuts for the final state
    par_S_profile{iterTau} = EntropyProfile(Lambda);

    % Correlations of the output field
    a=diag((1:dimB-1),1)/sqrt(dt);
    par_auto_corr_1{iterTau} = AutoCorrelation1(Gamma,Lambda,a',a,length(Gamma)-ceil(tau/dt),length(Gamma)-ceil(tau/dt)-round(CorrelationTime/dt));
    par_auto_corr_2{iterTau} = AutoCorrelation2(Gamma,Lambda,a',a,a',a,length(Gamma)-ceil(tau/dt),length(Gamma)-ceil(tau/dt)-round(CorrelationTime/dt));

end

excitations{iterphi}=par_excitations;
Lambdaloop{iterphi}=par_Lambdaloop;
S_profile{iterphi}=par_S_profile;
auto_corr_1{iterphi}=par_auto_corr_1;
auto_corr_2{iterphi}=par_auto_corr_2;
Purity{iterphi}=par_Purity;
PhotonsOut{iterphi}=par_PhotonsOut;
PhotonsIn{iterphi}=par_PhotonsIn;
PhotonNumberLoop{iterphi}=par_PhotonNumberLoop;

end
 
delete(gcp)

if saveresults==1
    save([datestr(datetime('now')) '.mat'])
end


if showfigures==1
%%
if length(phiInput)>1 && length(Tau)>1
    PurityScan=zeros(length(phiInput),length(Tau));
    ExcitationGround=zeros(length(phiInput),length(Tau));
    ExcitationLeft=zeros(length(phiInput),length(Tau));
    ExcitationRight=zeros(length(phiInput),length(Tau));
    PhotonEm=zeros(length(phiInput),length(Tau));
    PhotonInj=zeros(length(phiInput),length(Tau));
    EntropyLoop=zeros(length(phiInput),length(Tau));
    g2phi=zeros(length(phiInput),round(CorrelationTime/dt));
    g2tau=zeros(length(Tau),round(CorrelationTime/dt));
    
    for i=1:length(phiInput)
        for j=1:length(Tau)
            
            PurityScan(i,j)=Purity{i}{j}(end);
            ExcitationGround(i,j)=excitations{i}{j}(1,end);
            ExcitationRight(i,j)=excitations{i}{j}(2,end);
            ExcitationLeft(i,j)=excitations{i}{j}(3,end);
            PhotonEm(i,j)=real(PhotonsOut{i}{j}(end))/dt;
            PhotonInj(i,j)=real(PhotonsIn{i}{j}(end))/dt;
            EntropyLoop(i,j)=Lambdaloop{i}{j}{end}(1);
        end
        for j=1:round(CorrelationTime/dt)
            g2phi(i,j)=real(auto_corr_2{i}{end}(j)/(auto_corr_1{i}{end}(1))^2);
        end
    end
    for i=1:length(Tau)
        for j=1:round(CorrelationTime/dt)
            g2tau(i,j)=real(auto_corr_2{1}{i}(j)/(auto_corr_1{1}{i}(1))^2);
        end
    end
    
    figure
    subplot(3,3,1)
    [~, h2]=contourf(Tau,phiInput/pi,real(PurityScan),30);
    set(h2,'edgecolor','none')
    xlabel('\tau')
    ylabel('\phi')
    title('State purity')
    colormap('jet')
    colorbar
    subplot(3,3,2)
    [~, h2]=contourf(Tau,phiInput/pi,ExcitationGround,30);
    set(h2,'edgecolor','none')
    xlabel('\tau')
    ylabel('\phi')
    title('Ground state population')
    colorbar
    subplot(3,3,3)
    [~, h2]=contourf(Tau,phiInput/pi,ExcitationLeft,30);
    set(h2,'edgecolor','none')
    xlabel('\tau')
    ylabel('\phi')
    title('Left state population')
    colorbar
    subplot(3,3,4)
    [~, h2]=contourf(Tau,phiInput/pi,ExcitationRight,30);
    set(h2,'edgecolor','none')
    xlabel('\tau')
    ylabel('\phi')
    title('Right state population')
    colorbar
    subplot(3,3,5)
    [~, h2]=contourf(Tau,phiInput/pi,PhotonEm,30);
    set(h2,'edgecolor','none')
    xlabel('\tau')
    ylabel('\phi')
    title('Photons emmited per unit time')
    colorbar
    subplot(3,3,6)
    [~, h2]=contourf(Tau,phiInput/pi,PhotonInj,30);
    set(h2,'edgecolor','none')
    xlabel('\tau')
    ylabel('\phi')
    title('Photons injected in the loop per unit time')
    colorbar
    subplot(3,3,7)
    [~, h2]=contourf(Tau,phiInput/pi,EntropyLoop,30);
    set(h2,'edgecolor','none')
    xlabel('\tau')
    ylabel('\phi')
    title('Entropy loop-environment')
    colorbar
    subplot(3,3,8)
    [~, h2]=contourf((0:round(CorrelationTime/dt)-1)*dt,phiInput/pi,g2phi,30);
    set(h2,'edgecolor','none')
    xlabel('t')
    ylabel('\phi')
    title('g2(t)')
    colorbar
    subplot(3,3,9)
    [~, h2]=contourf((0:round(CorrelationTime/dt)-1)*dt,Tau,g2tau,30);
    set(h2,'edgecolor','none')
    xlabel('t')
    ylabel('\tau')
    title('g2(t)')
    colorbar
end
end