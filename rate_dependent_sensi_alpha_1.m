clear;

%% input variables
%phase numbers, the last one is the matrix
% load('slide4.mat','slideoutput');
% load('boundary4.mat','ttboundary');
%histogram(slideoutput(:,3),20)

%cracknumber=size(slideoutput,1);

%slideoutput(:,3)=30;

inputDS=[1*10^(-9) 5*10^(-9) 1*10^(-8) 5*10^(-8) 1*10^(-7)];
inputc=[0.01 0.02 0.04 0.08];
inputCK=[0.7 0.75 0.8 0.85 0.9];
inputalpha=[pi/18 pi/6 pi/4 pi/3];
sampling=20;
psampling=1;
ntheta=1;
phasen=psampling*sampling*ntheta+1;
inputtheta=[pi/18 pi/6 pi/4];


cstep(1)=1000;
cstep(2)=1000;

vstrain_c1=zeros(4,cstep(1)*2);


hstrain_c1=zeros(4,cstep(1)*2);


vstress_c1=zeros(4,cstep(1)*2);

for ai=1:4
    ai
    clearvars -except phasen inputalpha inputc inputCG inputCK...
        ai sampling  ntheta psampling inputtheta vstrain_c1 vstrain_c2 ...
        vstress_c1 vstress_c2 hstrain_c1 hstrain_c2 
%% parameters needed to be calibrated: check DS, RATE
uniformangel=linspace(0,360-360/sampling,sampling);
% chemial
% DS=inputDS(ai);  %mm^3/s
DS=2*10^(-10);  %mm^3/s
% mechanical
Km=17000;
Gm=10500;
nu=(3*Km-2*Gm)/2/(3*Km+Gm);
CK=0.9;
%CK=inputCK(ai);
 CG=0.6;
% CG=inputCG(ai);
% crack surface geometry
 alpha=inputalpha(ai);

 R13=1;
 R44=1;
 beta2=1.34;
 beta1=1.91;
% alpha=inputalpha(ai);

 c=0.002;
% c=inputc(ai);

theta_s=pi/6;
pr_s=0.5;
%%
% load experiment results
% load('bc8.mat');

% boundary=zeros(size(ttboundary,1),1);
% bi=0;
% for i = 1:size(ttboundary,1)
%     if ttboundary{i,1}>5
%         bi=bi+1;
% %        boundary(i)=ttboundary{i,1};
%     end
% end
%totalboundaries=sum(boundary);
% Em=25000; %young modulus mpa
% nu=0.32; % poisson ratio
% Km=Em/3/(1-2*nu);
% Gm=Em/2/(1+nu);

% invHC=inversensym2(HC);
sigma_o=[1 0 0; 0 1 0; 0 0 1]; % current global stress
sigma_devi=[0 0 0; 0 0 0; 0 0 6]; % global total stress
% rate=[0.1 0.01 0.001];
cstep=zeros(1,3);
rate=zeros(1,3);
cstep(1)=1000;
%cstep(2)=1000;
%cstep(3)=1022;
leisure=zeros(1,3);
leisure(1)=147;
leisure(2)=538;
 leisure=zeros(1,3);
cstep_ll=zeros(1,3);
tstep=0;
for i=1:2
    rate(i)=abs(sigma_devi(3,3)/cstep(i));
    cstep_ll(i)=cstep(i)*2+leisure(i);
    tstep=cstep_ll(i)+tstep;
end

inc=0.3;

ina=inc*1.01;
% crack direction


%thetao=zeros(size(slideoutput,1),1);
% phi=2*pi*rand(phasen,1);
% psi=2*pi*rand(phasen,1);

phio=pi/180*uniformangel;
psio=pi/180*uniformangel;


theta=zeros(phasen,1);
psi=zeros(phasen,1);
phi=zeros(phasen,1);

% for i=1:size(slideoutput,1)
%     thetao(i)=slideoutput(i,3)/180*pi;
% end 
% thetao=sort(thetao);
% 
% % Define the block parameter.  Average in a 100 row by 1 column wide window.
% blockSize = [round(size(slideoutput,1)/ntheta), 1];
% % Block process the image to replace every element in the 
% % 100 element wide block by the mean of the pixels in the block.
% % First, define the averaging function for use by blockproc().
% meanFilterFunction = @(theBlockStructure) mean2(theBlockStructure.data(:));
% % Now do the actual averaging (block average down to smaller size array).
% thetao = blockproc(thetao, blockSize, meanFilterFunction);


for i=1:psampling
    for j=1:sampling
        for k=1:ntheta
            theta((i-1)*sampling*ntheta+(j-1)*ntheta+k)=theta_s;
            psi((i-1)*sampling*ntheta+(j-1)*ntheta+k)=psio(j);
            phi((i-1)*sampling*ntheta+(j-1)*ntheta+k)=0;
        end
    end
end
psi(phasen)=0;
theta(phasen)=0;
phi(phasen)=0;

Q=cell(1,phasen);

for i=1:phasen
    Q(i)={transmatrixo(psi(i), theta(i), phi(i))};
end 



n0=zeros(phasen-1,3);
nx=zeros(phasen-1,3);
n1=zeros(phasen-1,3);
n2=zeros(phasen-1,3);

for i=1:phasen-1
    n0(i,:)=Q{i}'*[0 0 1]';
    nx(i,:)=Q{i}'*[1 0 0]';
    n1(i,:)=-sin(alpha)*nx(i,:)+cos(alpha)*n0(i,:);
    n2(i,:)=sin(alpha)*nx(i,:)+cos(alpha)*n0(i,:);
end 

%b=c/sin(alpha);



% chemical constants

C=6.46*10^(-6); % mol/mm^3
Omega=2.7*10^(4); % mm^3/mol
R=8.314; % J/K/mol
T=297; % K
V_d_constant=8*DS*C*Omega^2/R/T*sin(alpha)^3/c^2/1000;

% V_d_constant=0;
% inclusions shapes
% a=2:1:phasen+1; % oblate's radius
% c=0.5:1:phasen-0.5; % oblate's apurture
% a=[1 1 1 1 1 1 1 1 1 999999]
% c=2/3*a;
% c=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 999998.9999999]


% pstrain
pstrain=cell(1,phasen);
pstress=cell(1,phasen);
mstrain=cell(1,phasen);
lstrain=cell(1,phasen);

repstress=cell(1,tstep);
leta=cell(1,phasen);
letastrain=cell(1,phasen);
lvstrain=zeros(1,phasen);
terminate=0;
%% boundary conditions
% external stress

tstress=sigma_o;
vstress=zeros(tstep,3);
vstress(1,1)=1;
vstress(1,2)=1;
vstress(1,3)=1;
% external strain
% initial_strain=doubledotft(inversensym2(matrixc(Km,Gm)),sigma_o);
tstrain=cell(1,tstep);
vstrain=zeros(tstep,2);
hstrain=zeros(tstep,2);
% vstrain(1,1)=initial_strain(3,3);
% hstrain(1,1)=initial_strain(1,1);
for i=1:tstep
    tstrain(i)={zeros(3,3)};  
end  

% stiffness
pC=cell(1,phasen);

% volume ratio
pr=zeros(phasen,1);
for i=1:phasen-1
   pr(i)=pr_s/sampling;
end 
temp1=sum(pr(:));
pr(phasen)=1-temp1;

% % for k effect bc
% keffect=0.3;
% tstress=[-keffect*depth*0.0275 0 0; 0 -keffect*depth*0.0275 0; 0 0 -depth*0.0275];

%%
sigma_n1=zeros(phasen-1,1);
sigma_n2=zeros(phasen-1,1);

for i=1:phasen
    leta(i)={zeros(3,3)};  
    letastrain(i)={zeros(3,3)};
end 

% volumetric strain
tstrainv=zeros(1,tstep);

 % eigen strain
eta=zeros(1,phasen);
reta=zeros(tstep,1);
re_sigma_n1=zeros(tstep,1);
re_sigma_n2=zeros(tstep,1);

%% constants
% local stiffness

    pC{phasen}=matrixc(Km,Gm);


for i=1:phasen-1
%     pC{i}=transmatrix(psi(i), theta(i), phi(i),pC{i});
%     pC{i}(3,3)=0;
%     pC{i}=transmatrix(-psi(i), -theta(i), -phi(i),pC{i});
    pC{i}=matrixc(CK*Km,CG*Gm);
%    pC{i}=matrixc_traniso(CK*Km,CG*Gm,R13,R44);    
end


% P tensor
poP=cell(phasen,1); 
pP=cell(phasen,1); 
for i=1:phasen
    poP{i}=Ptensor(ina,inc,Km,Gm);  
    pP{i}=transmatrix(psi(i), theta(i), phi(i),poP{i});
end
% pP{10}=pP{1};
% infinite concentration tensor
pAo=cell(phasen,1); 
for i=1:phasen

    pAo{i}=inversensym2(symidendityf+doubledotff(pP{i},(pC{i}-pC{phasen})));
end



% average of infinite conentration tensor
sAo(1:3,1:3,1:3,1:3) = 0;
for i=1:phasen
    sAo=sAo+pr(i)*pAo{i};
end
% concentration tensor
pA=cell(phasen,1); 
for i=1:phasen
    pA{i}=doubledotff(pAo{i},inversensym2(sAo));
end

% homogenizaed stiffness tensor, Biot coefficient 
HC=zeros(3,3,3,3);
for i=1:phasen
    HC=pr(i)*doubledotff(pC{i},pA{i})+HC;
end

inversensym2(pA{1});

% biot modulus
AAP(1:3,1:3,1:3,1:3) = 0;
for i=1:phasen
    AAP=AAP+pr(i)*doubledotff(pAo{i},pP{i});
end
ACAP(1:3,1:3,1:3,1:3) = 0;
for i=1:phasen
    ACAP=ACAP+pr(i)*(doubledotff(doubledotff((HC-pC{i}),pAo{i}),pP{i}));
end
iN=cell(phasen,phasen); % inverse of biot modulus
for i=1:phasen
    for j=1:phasen
if i~=j 
iN{i,j}=pr(i)*(-doubledotff(doubledotff((pr(j)*pA{i}),pAo{j}),pP{j})+...
    doubledotff(doubledotff((doubledotff(pA{i},AAP)-doubledotff(pAo{i},pP{i})),inversensym2(ACAP)),...
     pr(j)*(transposefour(symidendityf-pA{j})+doubledotff(doubledotff((HC-pC{j}),pAo{j}),pP{j})))); 
else
iN{i,i}=pr(i)*(doubledotff(doubledotff((symidendityf-pr(i)*pA{i}),pAo{i}),pP{i})+...
    doubledotff(doubledotff((doubledotff(pA{i},AAP)-doubledotff(pAo{i},pP{i})),inversensym2(ACAP)),...
     pr(i)*(transposefour(symidendityf-pA{i})+doubledotff(doubledotff((HC-pC{i}),pAo{i}),pP{i})))); 
end
    end
end

% D tensor
Dt=cell(phasen,phasen);
 for i=1:phasen
    for j=1:phasen
     Dt{i,j}=doubledotff((iN{i,j}/pr(i)),pC{j});    
    end
 end
 
% parpool(4)
 
%% loading
tt=0;
for cycles=1:1
    ltt=0;
for ct=1:cstep_ll(cycles)
tt=tt+1;
ltt=ltt+1;




% weathering rate


for i=1:phasen-1
%    V_d=V_d_constant*(sigma_n1(i)-sigma_n2(i))/b/b;
% crack shear

    eta(i)=-V_d_constant*(sigma_n1(i)-sigma_n2(i))/inc;
%    
    if ct == 1
        eta(i)=0;
    end
    
    leta{i}(1,3)=leta{i}(1,3)+eta(i);
    leta{i}(3,1)=leta{i}(1,3);
    %
%     leta{i}(3,3)=leta{i}(3,3)-eta(i); 
    % chemical control
%    c(i)=c(i)*(1+eta(i));
    %
end

tletas=zeros(3,3);
for i=1:phasen-1
    letastrain(i)={Q{i}'*leta{i}*Q{i}};  
    tletas=tletas+letastrain{i};
end 



 % rev strain
 tempt2=zeros(3,3);
 for i=1:phasen
 tempt2=tempt2+pr(i)*doubledottf(doubledotft(pC{i},letastrain{i}),pA{i});
 end
 
 tstrain(tt)={doubledotft(inversensym2(HC),(tstress+tempt2))};
 
 % phase strain
  for j=1:phasen
  tempt1=zeros(3,3);
 for i=1:phasen
 tempt1=tempt1+doubledotft(Dt{j,i},letastrain{i});
 end
 pstrain{j}=doubledotft(pA{j},tstrain{tt})+tempt1;

 lstrain(j)={Q{j}*pstrain{j}*(Q{j})'};
 end
 
 % phase stress
  for j=1:phasen-1
 pstress{j}=doubledotft(pC{j},pstrain{j}-letastrain{j});
 sigma_n1(j)=heaviside(n1(j,:)*pstress{j}*n1(j,:)')*n1(j,:)*pstress{j}*n1(j,:)';
  sigma_n2(j)=heaviside(n2(j,:)*pstress{j}*n2(j,:)')*n2(j,:)*pstress{j}*n2(j,:)';
  end
  
 
repstress{tt}= pstress{1};
reta(ltt,cycles)=leta{1}(1,3);  
re_sigma_n1(ltt,cycles)= sigma_n1(1);  
re_sigma_n2(ltt,cycles)= sigma_n2(1);  
vstress(ltt,cycles)=tstress(3,3);
vstrain(ltt,cycles)=tstrain{tt}(3,3);
hstrain(ltt,cycles)=tstrain{tt}(1,1);

if ct<=cstep(cycles)
    tstress(3,3)=tstress(3,3)+rate(cycles);
elseif ct<=cstep(cycles)*2
    tstress(3,3)=tstress(3,3)-rate(cycles); 
end



    
end
end

vstrain_c1(ai,:)=vstrain(1:cstep_ll(1),1)-vstrain(1,1);

vstress_c1(ai,:)=vstress(1:cstep_ll(1),1)-vstress(1,1);

hstrain_c1(ai,:)=hstrain(1:cstep_ll(1),1)-hstrain(1,1);

end 

%% figure 



 figure
cruves{1}=plot(vstrain_c1(1,:),vstress_c1(1,:),'-k','LineWidth',1.5);
hold on

cruves{2}=plot(vstrain_c1(2,:),vstress_c1(2,:),'--k','LineWidth',1.5);
hold on

cruves{3}=plot(vstrain_c1(3,:),vstress_c1(3,:),'-.k','LineWidth',1.5);
hold on

cruves{4}=plot(vstrain_c1(4,:),vstress_c1(4,:),':k','LineWidth',1.5);
hold on

set(gca,'FontSize',16);
set(gcf, 'position', [200 200 500 400]);
ylabel('Differential stress (MPa)','fontsize',18);
xlabel('\epsilon_a','fontsize',18);
legend([cruves{1},cruves{2},cruves{3},cruves{4}],...
    '\alpha = 10^o','\alpha = 30^o','\alpha = 45^o',...
    '\alpha = 60^o','Location','southeast')

 figure
cruves{1}=plot(hstrain_c1(1,:),vstress_c1(1,:),'-k','LineWidth',1.5);
hold on
cruves{2}=plot(hstrain_c1(2,:),vstress_c1(2,:),'--k','LineWidth',1.5);
hold on
cruves{3}=plot(hstrain_c1(3,:),vstress_c1(3,:),'-.k','LineWidth',1.5);
hold on
cruves{4}=plot(hstrain_c1(4,:),vstress_c1(4,:),':k','LineWidth',1.5);
hold on

set(gca,'FontSize',16);
set(gcf, 'position', [200 200 500 400]);
ylabel('Differential stress (MPa)','fontsize',16);
xlabel('\epsilon_r','fontsize',16);
legend([cruves{1},cruves{2},cruves{3},cruves{4}],...
    '\alpha = 10^o','\alpha = 30^o','\alpha = 45^o',...
    '\alpha = 60^o','Location','southwest')

%% work

incre=size(hstrain_c1,2)-1;
rtw=zeros(4,incre);
tw_stress=zeros(4,incre);

for i=1:4
    tw=0;
    for j=1:incre
        tw_stress(i,j)=(vstress_c1(i,j)+vstress_c1(i,j+1))/2;
        dtw=(tw_stress(i,j)+1)*(vstrain_c1(i,j+1)-vstrain_c1(i,j))+2*(hstrain_c1(i,j+1)-hstrain_c1(i,j));
        tw=tw+dtw;
        rtw(i,j)=tw;
    end
    
end

 figure
cruves{1}=plot(tw_stress(1,:),rtw(1,:)/1000,'-k','LineWidth',1.5);
hold on
cruves{2}=plot(tw_stress(2,:),rtw(2,:)/1000,'--k','LineWidth',1.5);
hold on
cruves{3}=plot(tw_stress(3,:),rtw(3,:)/1000,'-.k','LineWidth',1.5);
hold on
cruves{4}=plot(tw_stress(4,:),rtw(4,:)/1000,':k','LineWidth',1.5);
hold on

set(gca,'FontSize',16);
set(gcf, 'position', [200 200 500 400]);
xlabel('Differential stress (MPa)','fontsize',16);
ylabel('External unit work (J/mm^3)','fontsize',16);
legend([cruves{1},cruves{2},cruves{3},cruves{4}],...
'\alpha = 10^o','\alpha = 30^o','\alpha = 45^o',...
    '\alpha = 60^o','Location','northwest')


% 
% figure
% plot(cftstrainv);
% set(gca,'FontSize',16);
% set(gcf, 'position', [200 200 400 300]);
% ylabel('cstrainv','fontsize',16);
% xlabel('time(years)','fontsize',16);
% 
% 
% figure
% plot(fomega)
% set(gca,'FontSize',16);
% set(gcf, 'position', [200 200 400 300]);
% xlabel('time(years)','fontsize',16);
% ylabel('omega','fontsize',16);
% 
% figure
% plot(fc1)
% set(gca,'FontSize',16);
% set(gcf, 'position', [200 200 400 300]);
% xlabel('time(years)','fontsize',16);
% ylabel('c1(mm)','fontsize',16);
% 
% figure
% plot(fc2)
% set(gca,'FontSize',16);
% set(gcf, 'position', [200 200 400 300]);
% xlabel('time(years)','fontsize',16);
% ylabel('c2(mm)','fontsize',16);
% 
% figure
% plot(fc3)
% set(gca,'FontSize',16);
% set(gcf, 'position', [200 200 400 300]);
% xlabel('time(years)','fontsize',16);
% ylabel('c3(mm)','fontsize',16);
% figure
% plot(fa1)
% set(gca,'FontSize',16);
% set(gcf, 'position', [200 200 400 300]);
% xlabel('time(years)','fontsize',16);
% ylabel('a1(mm)','fontsize',16);
% figure
% plot(fa2)
% set(gca,'FontSize',16);
% set(gcf, 'position', [200 200 400 300]);
% xlabel('time(years)','fontsize',16);
% ylabel('a2(mm)','fontsize',16);
% figure
% plot(fa3)
% set(gca,'FontSize',16);
% set(gcf, 'position', [200 200 400 300]);
% xlabel('time(years)','fontsize',16);
% ylabel('a3(mm)','fontsize',16);