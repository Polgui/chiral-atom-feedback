clear all
n=4;


sm=sparse([0,0;1,0]);
sx=sparse([0,1;1,0]);
sy=sparse([0,-1i;1i,0]);
sz=sparse([1,0;0,-1]);
id=speye(2);
SWAP=sparse([1,0,0,0;0,0,1,0;0,1,0,0;0,0,0,1]);

SX=cell(n,1);
SY=cell(n,1);
SZ=cell(n,1);
SM=cell(n,1);
for(k=1:n)
    SX{k}=KroneckerProduct(speye(2^(k-1)),sx,speye(2^(n-k)));
    SY{k}=KroneckerProduct(speye(2^(k-1)),sy,speye(2^(n-k)));
    SZ{k}=KroneckerProduct(speye(2^(k-1)),sz,speye(2^(n-k)));
    SM{k}=KroneckerProduct(speye(2^(k-1)),sm,speye(2^(n-k)));
end

Jx=spalloc(2^n,2^n,1);
Jy=spalloc(2^n,2^n,1);
Jz=spalloc(2^n,2^n,1);
Jm=spalloc(2^n,2^n,1);
for k=1:n
    Jx=Jx+SX{k};
    Jy=Jy+SY{k};
    Jz=Jz+SZ{k};
    Jm=Jm+SM{k};
end

Omegapattern=0*ones(1,n);
%Omegapattern=[1,0,-1,0];
deltapattern=0*ones(1,n);
%deltapattern=[0,1/8,-1/8,0];

Hsys=spalloc(2^n,2^n,1);
for k=1:n
    Hsys=Hsys+Omegapattern(k)*SX{k}-deltapattern(k)*SZ{k};
end

Ldrive=-1i*SOpre(Hsys)+1i*SOpost(Hsys);

% Hcasc=spalloc(2^n,2^n,1);
% for k=1:n
%     for l=k+1:n
%         Hcasc=Hcasc-1i/2*(SM{l}'*SM{k}-SM{l}*SM{k}');
%     end
% end

channel{1}=[1,2];
channel{2}=[2,3];
channel{3}=[3,4];
channel{4}=[4,1];
channel{5}=[2,1];
channel{6}=[3,2];
channel{7}=[4,3];
channel{8}=[1,4];

gammachannel{1}=1;
gammachannel{2}=1;
gammachannel{3}=1;
gammachannel{4}=1;
gammachannel{5}=1;
gammachannel{6}=1;
gammachannel{7}=1;
gammachannel{8}=1;

Nchannel=size(channel,2);
Hchannel=cell(1,Nchannel);
cchannel=cell(1,Nchannel);
Lchannel=cell(1,Nchannel);

Ltotal=spalloc(2^(2*n),2^(2*n),1);
for m=1:Nchannel
Hchannel{m}=spalloc(2^n,2^n,1);
cchannel{m}=spalloc(2^n,2^n,1);
    for k=1:size(channel{m},2);
        site1=channel{m}(k);
        for l=k+1:size(channel{m},2)
            site2=channel{m}(l);
            Hchannel{m}=Hchannel{m}-1i/2*(SM{site2}'*SM{site1}-SM{site2}*SM{site1}');
        end
        cchannel{m}=cchannel{m}+SM{site1};
    end
D=SOpre(cchannel{m})*SOpost(cchannel{m}')-1/2*SOpre(cchannel{m}'*cchannel{m})-1/2*SOpost(cchannel{m}'*cchannel{m});
Lchannel{m}=gammachannel{m}*(-1i.*SOpre(Hchannel{m})+1i*SOpost(Hchannel{m})+D);
Ltotal=Ltotal+Lchannel{m};
end
Ltotal=Ltotal+Ldrive;

Groundstate=KroneckerProduct([0,1],[0,1],[0,1],[0,1])';
S=(SM{1}'-SM{2}'+SM{3}'-SM{4}')*Groundstate;

rhoS=S*S';
rhoSvec=rhoS(:);


%[EV,EW]=eigs(Ltotal,2,'sm');
[EV,EW]=eig(full(Ltotal));
[EW,ind]=sort(abs(diag(EW)),'ascend');
EV=EV(:,ind);
liouvill_Gap=EW(2)-EW(1)

rhoSSvec=EV(:,1);
%rhoSS=vec2mat(rhoSSvec);
rhoSS=reshape(rhoSSvec,2^n,2^n);
rhoSS=rhoSS./trace(rhoSS);rhoSS=(rhoSS+rhoSS')./2;
Purity=trace(rhoSS*rhoSS)

% [EV,EW]=eig(full(rhoSS));
% [EW,ind]=sort(diag(EW));
% psi=EV(:,end);psi=psi./norm(psi);
% % rhoSS12=TrX(rhoSS,[3,4],[2,2,2,2]);S12=Entropy(rhoSS12)
% % rhoSS13=TrX(rhoSS,[2,4],[2,2,2,2]);S13=Entropy(rhoSS13)
% % rhoSS14=TrX(rhoSS,[2,3],[2,2,2,2]);S14=Entropy(rhoSS14)
% % rhoSS23=TrX(rhoSS,[1,4],[2,2,2,2]);S23=Entropy(rhoSS23)
% % rhoSS24=TrX(rhoSS,[1,3],[2,2,2,2]);S24=Entropy(rhoSS24)
% % rhoSS34=TrX(rhoSS,[1,2],[2,2,2,2]);S34=Entropy(rhoSS34)
% 
% 


SvN=Entropy(rhoSS)
%FQG_SS=Quantum_fisher_info(rhoSS,Generator)
%FC_SS=Fisher_info_mixed(full(rhoSS),full(Generator),full(Jz))

dim=2*ones(1,n);
sys=[1:n];

for(k=1:n-1)
    for(l=k+1:n)
        sys=[1:n];
        sys(l)=[];
        sys(k)=[];
    rhoSSkl=TrX(rhoSS,sys,dim);
    S_pair(k,l)=Entropy(rhoSSkl);
    end
end
S_pair
%%%%%%%%%%%%%%%%% TIME EVOLUTION %%%%%%%%%%%%%%%%%%%%%%

rho0=eye(2^n)./2^n;
%rho0=zeros(2^n);rho0(end,end)=1;
rho0vec=mat2vec(rho0);

Nt=80;
tspan=linspace(0,20,Nt);dt=tspan(2)-tspan(1);
rhovect(:,1)=rho0vec;
%SOURCE=zeros(Nt,n);
for k=1:Nt
    rho=vec2mat(rhovect(:,k));rho=(rho+rho')./2;rho=rho./trace(rho);
    Purity_t(k)=trace(rho*rho);
%     Entropy_t(k)=Entropy(rho);
%     for(l=1:n-1)
%         sys1=sys;
%         sys1(l:l+1)=[];
%         rhoreduced=TrX(rho,[sys1],dim);
%         
%         S(l,k)=Entropy(rhoreduced);
%     end
    
%     for l=1:n
%         SOURCE(k,l)=trace(rho*SM{l});
%     end
%  
%     Eleft(k,1)=0;
%     Eright(k,n+1)=0;
%     for l=2:n+1
%         Eleft(k,l)=1i.*sqrt(GammaL)*sum(SOURCE(k,1:l-1));
%         Eright(k,l-1)=1i.*sqrt(GammaR)*sum(SOURCE(k,l-1:n));
%     end
%     
    
    rhovect(:,k+1)=expv(dt,Ltotal,rhovect(:,k));
end

figure(1)
plot(tspan,Purity_t)
% % figure(2)
% % plot(tspan,FQG,'b','Linewidth',1.1)
% % hold on
% % plot(tspan,FC,'r','Linewidth',1.1)
% % hold off
% % hold on
% %figure(3)
% %plot(tspan,C)
% %figure(4)
% %imagesc(real(Eright+Eleft).^2)
% % 
% % figF=figure(2)
% % set(gcf, 'PaperUnits', 'centimeters');
% % set(gcf, 'PaperSize', [9 6]);
% % set(gcf, 'PaperPosition', [0 0 9 6]);
% % print('-dpdf','-r300',[figF,'figF.pdf'])
