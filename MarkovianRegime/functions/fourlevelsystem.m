clear all
n=4;


smminus=sparse([0,0,0,0;0,0,0,0;0,0,0,0;1,0,0,0]);
smplus=sparse([0,0,0,0;0,0,0,0;0,0,0,0;0,0,1,0]);
jm=sparse([0,1,0,0;0,0,1,0;0,0,0,0;0,0,0,0]);
sx=sparse([0,1;1,0]);
sy=sparse([0,-1i;1i,0]);
sz=sparse([1,0;0,-1]);
id=speye(4);


SMminus=cell(n,1);
SMplus=cell(n,1);
for(k=1:n)
    SMminus{k}=KroneckerProduct(speye(4^(k-1)),smminus,speye(4^(n-k)));
    SMplus{k}=KroneckerProduct(speye(4^(k-1)),smplus,speye(4^(n-k)));
    JM{k}=KroneckerProduct(speye(4^(k-1)),jm,speye(4^(n-k)));
end
B=100*ones(1,n);
Omegapattern=0*ones(1,n);
%Omegapattern=[1,0,-1,0];
deltapattern=0*ones(1,n);
%deltapattern=[0,1/8,-1/8,0];

Hsys=spalloc(4^n,4^n,1);
for k=1:n
    Hsys=Hsys+B(k)*(JM{k}+JM{k}');
end

Ldrive=-1i*SOpre(Hsys)+1i*SOpost(Hsys);

% Hcasc=spalloc(2^n,2^n,1);
% for k=1:n
%     for l=k+1:n
%         Hcasc=Hcasc-1i/2*(SM{l}'*SM{k}-SM{l}*SM{k}');
%     end
% end

channelminus{1}=[1,2];
channelminus{2}=[2,3];
channelminus{3}=[3,4];
channelminus{4}=[4,1];
channelplus{1}=[2,1];
channelplus{2}=[3,2];
channelplus{3}=[4,3];
channelplus{4}=[1,4];

gammachannel{1}=1;
gammachannel{2}=1;
gammachannel{3}=1;
gammachannel{4}=1;


Nchannel=4;
Hchannel=cell(1,Nchannel);
cchannel=cell(1,Nchannel);
Lchannel=cell(1,Nchannel);

Ltotal=spalloc(4^(2*n),2^(4*n),1);
for m=1:Nchannel
Hchannelminus{m}=spalloc(4^n,4^n,1);
cchannelminus{m}=spalloc(4^n,4^n,1);
    for k=1:size(channelminus{m},2);
        site1=channelminus{m}(k);
        for l=k+1:size(channelminus{m},2)
            site2=channelminus{m}(l);
            Hchannelminus{m}=Hchannelminus{m}-1i/2*(SMminus{site2}'*SMminus{site1}-SMminus{site2}*SMminus{site1}');
        end
        cchannelminus{m}=cchannelminus{m}+SMminus{site1};
    end
D=SOpre(cchannelminus{m})*SOpost(cchannelminus{m}')-1/2*SOpre(cchannelminus{m}'*cchannelminus{m})-1/2*SOpost(cchannelminus{m}'*cchannelminus{m});
Lchannelminus{m}=gammachannel{m}*(-1i.*SOpre(Hchannelminus{m})+1i*SOpost(Hchannelminus{m})+D);

Hchannelplus{m}=spalloc(4^n,4^n,1);
cchannelplus{m}=spalloc(4^n,4^n,1);
    for k=1:size(channelplus{m},2);
        site1=channelplus{m}(k);
        for l=k+1:size(channelplus{m},2)
            site2=channelplus{m}(l);
            Hchannelplus{m}=Hchannelplus{m}-1i/2*(SMplus{site2}'*SMplus{site1}-SMplus{site2}*SMplus{site1}');
        end
        cchannelplus{m}=cchannelplus{m}+SMplus{site1};
    end
D=SOpre(cchannelplus{m})*SOpost(cchannelplus{m}')-1/2*SOpre(cchannelplus{m}'*cchannelplus{m})-1/2*SOpost(cchannelplus{m}'*cchannelplus{m});
Lchannelplus{m}=gammachannel{m}*(-1i.*SOpre(Hchannelplus{m})+1i*SOpost(Hchannelplus{m})+D);

Ltotal=Ltotal+Lchannelminus{m}+Lchannelplus{m};
end
Ltotal=Ltotal+Ldrive;

Groundstate=KroneckerProduct([0,0,0,1],[0,0,0,1],[0,0,0,1],[0,0,0,1])';
S=(SMminus{1}'-SMminus{2}'+SMminus{3}'-SMminus{4}'+SMplus{1}'-SMplus{2}'+SMplus{3}'-SMplus{4}')*Groundstate;

rhoS=S*S';
rhoSvec=rhoS(:);


%[EV,EW]=eigs(Ltotal,2,'sm');
% [EV,EW]=eig(full(Ltotal));
% [EW,ind]=sort(abs(diag(EW)),'ascend');
% EV=EV(:,ind);
% liouvill_Gap=EW(2)-EW(1)
% 
% rhoSSvec=EV(:,1);
% %rhoSS=vec2mat(rhoSSvec);
% rhoSS=reshape(rhoSSvec,2^n,2^n);
% rhoSS=rhoSS./trace(rhoSS);rhoSS=(rhoSS+rhoSS')./2;
% Purity=trace(rhoSS*rhoSS)
% 
% % [EV,EW]=eig(full(rhoSS));
% % [EW,ind]=sort(diag(EW));
% % psi=EV(:,end);psi=psi./norm(psi);
% % % rhoSS12=TrX(rhoSS,[3,4],[2,2,2,2]);S12=Entropy(rhoSS12)
% % % rhoSS13=TrX(rhoSS,[2,4],[2,2,2,2]);S13=Entropy(rhoSS13)
% % % rhoSS14=TrX(rhoSS,[2,3],[2,2,2,2]);S14=Entropy(rhoSS14)
% % % rhoSS23=TrX(rhoSS,[1,4],[2,2,2,2]);S23=Entropy(rhoSS23)
% % % rhoSS24=TrX(rhoSS,[1,3],[2,2,2,2]);S24=Entropy(rhoSS24)
% % % rhoSS34=TrX(rhoSS,[1,2],[2,2,2,2]);S34=Entropy(rhoSS34)
% % 
% % 
% 
% 
% SvN=Entropy(rhoSS)
% %FQG_SS=Quantum_fisher_info(rhoSS,Generator)
% %FC_SS=Fisher_info_mixed(full(rhoSS),full(Generator),full(Jz))
% 
% dim=2*ones(1,n);
% sys=[1:n];
% 
% for(k=1:n-1)
%     for(l=k+1:n)
%         sys=[1:n];
%         sys(l)=[];
%         sys(k)=[];
%     rhoSSkl=TrX(rhoSS,sys,dim);
%     S_pair(k,l)=Entropy(rhoSSkl);
%     end
% end
% S_pair
%%%%%%%%%%%%%%%%% TIME EVOLUTION %%%%%%%%%%%%%%%%%%%%%%

rho0=eye(4^n)./2^n;
%rho0=zeros(2^n);rho0(end,end)=1;
rho0vec=mat2vec(rho0);

Nt=20;
tspan=linspace(0,100,Nt);dt=tspan(2)-tspan(1);
rhovect(:,1)=rho0vec;
%SOURCE=zeros(Nt,n);
for k=1:Nt
    rho=vec2mat(rhovect(:,k));rho=(rho+rho')./2;rho=rho./trace(rho);
    Purity_t(k)=trace(rho*rho);
    rhovect(:,k+1)=expv(dt,Ltotal,rhovect(:,k));
end
% 
figure(1)
plot(tspan,Purity_t)
% % % figure(2)
% % % plot(tspan,FQG,'b','Linewidth',1.1)
% % % hold on
% % % plot(tspan,FC,'r','Linewidth',1.1)
% % % hold off
% % % hold on
% % %figure(3)
% % %plot(tspan,C)
% % %figure(4)
% % %imagesc(real(Eright+Eleft).^2)
% % % 
% % % figF=figure(2)
% % % set(gcf, 'PaperUnits', 'centimeters');
% % % set(gcf, 'PaperSize', [9 6]);
% % % set(gcf, 'PaperPosition', [0 0 9 6]);
% % % print('-dpdf','-r300',[figF,'figF.pdf'])
