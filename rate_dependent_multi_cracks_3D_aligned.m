clear;

%% input variables
%phase numbers, the last one is the matrix
load('slide4.mat','slideoutput');
load('boundary4.mat','ttboundary');
%histogram(slideoutput(:,3),20)

cracknumber=size(slideoutput,1);

%slideoutput(:,3)=30;

inputDS=[1*10^(-9) 5*10^(-9) 1*10^(-8) 5*10^(-8) 1*10^(-7)];
inputCG=[0.7 0.75 0.8 0.85 0.9];
inputCK=[0.7 0.75 0.8 0.85 0.9];
inputalpha=[pi/18 pi/9 pi/6 pi/4 pi/3];
sampling=10;

phasen=sampling*sampling*sampling+1;

for ai=1:1
    ai
    clearvars -except phasen slideoutput inputalpha inputDS inputCG inputCK ttboundary ai sampling cracknumber
%% parameters needed to be calibrated: check DS, RATE
uniformangel=linspace(0,360-360/sampling,sampling);
% chemial
% DS=inputDS(ai);  %mm^3/s
DS=1*10^(-9);  %mm^3/s
% mechanical
Km=30000;
Gm=14400;
nu=(3*Km-2*Gm)/2/(3*Km+Gm);
CK=0.9;
%CK=inputCK(ai);
 CG=0.8;
% CG=inputCG(ai);
% crack surface geometry
 alpha=pi/18;

% alpha=inputalpha(ai);

 c=0.02;
% c=inputc(ai);
%%
% load experiment results
load('bc8.mat');

% boundary=zeros(size(ttboundary,1),1);
bi=0;
for i = 1:size(ttboundary,1)
    if ttboundary{i,1}>5
        bi=bi+1;
%        boundary(i)=ttboundary{i,1};
    end
end
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
cstep(1)=10;
cstep(2)=100;
cstep(3)=1022;
leisure=zeros(1,3);
leisure(1)=147;
leisure(2)=538;
% leisure=zeros(1,3);
cstep_ll=zeros(1,3);
tstep=0;
for i=1:3
    rate(i)=abs(sigma_devi(3,3)/cstep(i));
    cstep_ll(i)=cstep(i)*2+leisure(i);
    tstep=cstep_ll(i)+tstep;
end

inc=0.3;

ina=inc*1.01;
% crack direction


thetao=zeros(size(slideoutput,1),1);
% phi=2*pi*rand(phasen,1);
% psi=2*pi*rand(phasen,1);

phio=pi/180*uniformangel;
psio=pi/180*uniformangel;


theta=zeros(size(phasen,1),1);
psi=zeros(size(phasen,1),1);
phi=zeros(size(phasen,1),1);

for i=1:size(slideoutput,1)
    thetao(i)=slideoutput(i,3)/180*pi;
end 
thetao=sort(thetao);

% Define the block parameter.  Average in a 100 row by 1 column wide window.
blockSize = [round(size(slideoutput,1)/sampling), 1];
% Block process the image to replace every element in the 
% 100 element wide block by the mean of the pixels in the block.
% First, define the averaging function for use by blockproc().
meanFilterFunction = @(theBlockStructure) mean2(theBlockStructure.data(:));
% Now do the actual averaging (block average down to smaller size array).
thetao = blockproc(thetao, blockSize, meanFilterFunction);


for i=1:sampling
    for j=1:sampling
        for k=1:sampling
            theta((i-1)*sampling*sampling+(j-1)*sampling+k)=thetao(i);
            psi((i-1)*sampling*sampling+(j-1)*sampling+k)=psio(j);
            phi((i-1)*sampling*sampling+(j-1)*sampling+k)=phio(k);
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
V_d_constant=8*DS*C*Omega^2/R/T*sin(alpha)/c^2/1000;

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
vstrain=zeros(tstep,3);
hstrain=zeros(tstep,3);
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
   pr(i)=1/sampling/sampling/sampling*((size(slideoutput,1))/bi)^(0.5);
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
reta=zeros(tstep,3);
re_sigma_n1=zeros(tstep,3);
re_sigma_n2=zeros(tstep,3);

%% constants
% local stiffness

    pC{phasen}=matrixc(Km,Gm);


for i=1:phasen-1
%     pC{i}=transmatrix(psi(i), theta(i), phi(i),pC{i});
%     pC{i}(3,3)=0;
%     pC{i}=transmatrix(-psi(i), -theta(i), -phi(i),pC{i});
    pC{i}=matrixc(CK*Km,CG*Gm);
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
    HC=pr(1)*doubledotff(pC{1},pA{1})+pr(phasen)*doubledotff(pC{phasen},pA{phasen});


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
        if j==1 || j==phasen
if i~=j 
iN{i,j}=pr(i)*(-doubledotff(doubledotff((pr(j)*pA{i}),pAo{j}),pP{j})+...
    doubledotff(doubledotff((doubledotff(pA{i},AAP)-doubledotff(pAo{i},pP{i})),inversensym2(ACAP)),...
     pr(j)*(transposefour(symidendityf-pA{j})+doubledotff(doubledotff((HC-pC{j}),pAo{j}),pP{j})))); 
else
iN{i,i}=pr(i)*(doubledotff(doubledotff((symidendityf-pr(i)*pA{i}),pAo{i}),pP{i})+...
    doubledotff(doubledotff((doubledotff(pA{i},AAP)-doubledotff(pAo{i},pP{i})),inversensym2(ACAP)),...
     pr(i)*(transposefour(symidendityf-pA{i})+doubledotff(doubledotff((HC-pC{i}),pAo{i}),pP{i})))); 
end
        else
            iN{i,j}=zeros(3,3,3,3);
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
for cycles=1:3
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


for i=1:phasen-1
    letastrain(i)={Q{i}'*leta{i}*Q{i}};  

end 



 % rev strain
 tempt1=zeros(3,3);
 for i=1:phasen
 tempt1=tempt1+pr(i)*doubledottf(doubledotft(pC{i},letastrain{i}),pA{i});
 end
 
 tstrain(tt)={doubledotft(inversensym2(HC),(tstress+tempt1))};
 
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


if tt==cstep(cycles)


end
    
end
end

%% figure 
load('bc8.mat'); 


figure('Name','sn','NumberTitle','off')
cruves{1}=plot(vstrain(1:cstep_ll(1),1)-vstrain(1,1),re_sigma_n1(1:cstep_ll(1),1),'-k','LineWidth',1);
hold on
cruves{2}=plot(vstrain(1:cstep_ll(2),2)-vstrain(1,1),re_sigma_n1(1:cstep_ll(2),2),'-r','LineWidth',1);
hold on
cruves{3}=plot(vstrain(1:cstep_ll(3),3)-vstrain(1,1),re_sigma_n1(1:cstep_ll(3),3),'-b','LineWidth',1);
hold on
cruves{4}=plot(vstrain(1:cstep_ll(1),1)-vstrain(1,1),re_sigma_n2(1:cstep_ll(1),1),'--k','LineWidth',1);
hold on
cruves{5}=plot(vstrain(1:cstep_ll(2),2)-vstrain(1,1),re_sigma_n2(1:cstep_ll(2),2),'--r','LineWidth',1);
hold on
cruves{6}=plot(vstrain(1:cstep_ll(3),3)-vstrain(1,1),re_sigma_n2(1:cstep_ll(3),3),'--b','LineWidth',1);
hold on

set(gcf, 'position', [200 200 500 400]);
set(gca,'FontSize',18)
xlabel('\epsilon_a','fontsize',18)
ylabel('\sigma_n','fontsize',18)
legend([cruves{1},cruves{2},cruves{3},cruves{4},cruves{5},cruves{6}],...
    '\sigma_{n1}1','\sigma_{n1}2','\sigma_{n1}3','\sigma_{n2}1','\sigma_{n2}2','\sigma_{n2}3','Location','best')
savefig(['rate_d_sn_CK',num2str(ai),'.fig'])
% 
% 
 figure
cruves{1}=plot(vstrain(1:cstep_ll(1),1)-vstrain(1,1),reta(1:cstep_ll(1),1),'-k','LineWidth',1);
hold on
cruves{2}=plot(vstrain(1:cstep_ll(2),2)-vstrain(1,1),reta(1:cstep_ll(2),2),'-r','LineWidth',1);
hold on
cruves{3}=plot(vstrain(1:cstep_ll(3),3)-vstrain(1,1),reta(1:cstep_ll(3),3),'-b','LineWidth',1);
hold on
set(gca,'FontSize',16);
set(gcf, 'position', [200 200 500 400]);
ylabel('\epsilon_c','fontsize',18);
xlabel('\epsilon_a','fontsize',18);
legend([cruves{1},cruves{2},cruves{3}],...
    'cycle 1','cycle 2','cycle 3','Location','southeast')
savefig(['rate_d_se_CK',num2str(ai),'.fig'])


 figure
cruves{1}=plot(vstrain(1:cstep_ll(1),1)-vstrain(1,1),vstress(1:cstep_ll(1),1)-vstress(1,1),'-k','LineWidth',2);
hold on
cruves{2}=plot(vstrain(1:cstep_ll(2),2)+bc{2}(1,2)-bc{1}(1,2)-vstrain(1,2),vstress(1:cstep_ll(2),2)-vstress(1,2),'-r','LineWidth',2);
hold on
cruves{3}=plot(vstrain(1:cstep_ll(3),3)+bc{3}(1,2)-bc{1}(1,2)-vstrain(1,3),vstress(1:cstep_ll(3),3)-vstress(1,3),'-b','LineWidth',2);
hold on
cruves{4}=plot(bc{1}(:,2)-bc{1}(1,2),bc{1}(:,5),'--k','LineWidth',2);
cruves{5}=plot(bc{2}(:,2)-bc{1}(1,2),bc{2}(:,5),'--r','LineWidth',2);
cruves{6}=plot(bc{3}(:,2)-bc{1}(1,2),bc{3}(:,5),'--b','LineWidth',2);

set(gca,'FontSize',16);
set(gcf, 'position', [200 200 500 400]);
ylabel('\sigma_a','fontsize',18);
xlabel('\epsilon_a','fontsize',18);
legend([cruves{1},cruves{2},cruves{3},cruves{4},cruves{5},cruves{6}],...
    'cycle 1-s','cycle 2-s','cycle 3-s','cycle 1-e','cycle 2-e','cycle 3-e','Location','southeast')

savefig(['rate_d_ssv_CK',num2str(ai),'.fig'])

 figure
cruves{1}=plot(hstrain(1:cstep_ll(1),1)-hstrain(1,1),vstress(1:cstep_ll(1),1)-vstress(1,1),'-k','LineWidth',2);
hold on
cruves{2}=plot(hstrain(1:cstep_ll(2),2)+bc{2}(1,3)-bc{1}(1,3)-hstrain(1,2),vstress(1:cstep_ll(2),2)-vstress(1,2),'-r','LineWidth',2);
hold on
cruves{3}=plot(hstrain(1:cstep_ll(3),3)+bc{3}(1,3)-bc{1}(1,3)-hstrain(1,3),vstress(1:cstep_ll(3),3)-vstress(1,3),'-b','LineWidth',2);
hold on
cruves{4}=plot(bc{1}(:,3)-bc{1}(1,3),bc{1}(:,5),'--k','LineWidth',2);
cruves{5}=plot(bc{2}(:,3)-bc{1}(1,3),bc{2}(:,5),'--r','LineWidth',2);
cruves{6}=plot(bc{3}(:,3)-bc{1}(1,3),bc{3}(:,5),'--b','LineWidth',2);

set(gca,'FontSize',16);
set(gcf, 'position', [200 200 500 400]);
ylabel('\sigma_a','fontsize',16);
xlabel('\epsilon_l','fontsize',16);
legend([cruves{1},cruves{2},cruves{3},cruves{4},cruves{5},cruves{6}],...
    'C1 (Simulation)','C2 (Simulation)','C3 (Simulation)','C1 (Experiment)','C2 (Experiment)','C3 (Experiment)','Location','northeast')

savefig(['rate_d_ssh_CK',num2str(ai),'.fig'])
%% test

end 

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