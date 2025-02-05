function defoliate
% Code for solving thornley model converted thornley's input values 
addpath('/home/bop15bh')
close ALL
fa=[];
fb=[];
alun=[];
kcc1=[];
kcc2=[];knn1=[];knn2=[];
N1=[];
N2=[];
C1=[];
C2=[];
%% input parameters
fly=40; % Carbon uptake rate units: micromol*m^(-2)*s^(-1)
theta=(1).*10.^(6); % "theta" dry density of a leaf units: g*m^(-3)
sky=61; %vmax for N uptake rate
Yg=0.66; % conversion efficieny from substrate to plant material
LL=1./10; % "lambda" N:C ratio in the leaf
LR=1./20; % "lambda" N:C ratio in the root
den=120; %leaf thickness units: micrometre used for carbon uptake rate
rif=730; % root thickness micrometres 
thick=.24.*10^(6); % root tissue density units: g*m^-3
% Concentration constraints in nmolmg^(-1)
lcmin0=92.8;
lnmin0=0.1;
rcmin0=63;
rnmin0=7.54;
k1=1000; k2=1000;  % k1,k2 and v1,v2 control RGR 
V1=60;
V2=60;  
r1=(15).*10.^(-12);% leaf respiration rate kgmolg^-1s^-1
r2=(10).*10.^(-12);% root respiration rate kgmolg^-1s^-1
ld=10.^6; %leaf tissue density gm^-3
rd=0.2.*10.^6; % root tissue density gm^-3
R1=r1.*ld
R2=r2.*rd

knn=103;
kcc=200;
%% Converting to thornley units/dimensions 
%-------------------------------------------------------------------------
% definitions
x0=(1./(theta)).*((10.^3).*12) % "theta" conversion factor for plant material to volume converting input into thornley units: m^3*(kgmol)^-1
k=9.*10.^(-4); % involved in rate of use of substrate
rc=0.5.*10.^(3); % carbon transport resistence
rn=1.*10.^(3); % nitrogen transport resistence
Kn=sky.*thick.*(10.^(-9).*10.^(-3)).*1 % Nitrogen uptake rate converting micromol*kg^(-1)*s^(-1) to kg*m^-3 by 10^-9 then 10^-9 from micromol to kgmol. thornley units: kgmol*m^-3s^-1
Kc=(10.^(-3).*fly./(den))%.*1.669.*1%./4.82; % Carbon uptake rate=photosynthesis/(leaf thickness in metres) converted from mumol to kgmol
x0C=theta.*0.45;%Theta conversion for c concentration
x0NL=theta.*0.45.*LL;% theta conversion for n concentration
x0NR=theta.*0.45.*LR;
% initial conditions 
ileaf=0.01;
iroot=0.01;
j=400.*x0C.*10.^(-9);
a1=10.^(-5);burger=10.^(-5);
%burger=0.45;
lam=Kc./Kn
% converting leaf concentrations into thornley units
lcmin=lcmin0.*x0C.*10.^(-9)
lnmin=lnmin0.*x0NL.*10.^(-9)
% converting root concentrations into thornley units
rcmin=rcmin0.*x0C.*10.^(-9)
rnmin=rnmin0.*x0NR.*10.^(-9)

%tspan=[0:200:4665600]; % time span for 2years

    C={'r','b','g','k','y','c','m'}%
        % switching the feedbacks on 
             aa=1;bb=1;cc=1;d1=1;ee=1;f1=1; 
n=400;% sets soil N concentration
% Setting atmospheric Co2 concentrations
for i=1:2
         if i==1
            Aco2=350; 
         else
        Aco2=700;         
         end
         co2=0.7.*Aco2; % intercellular co2

for o=1:2
 %% Solving the model for the first 10 days 
    if o==1
        tspan=[0:100:864000]; %time span for 20 days 
        y1_0=ileaf./theta;
        y2_0=iroot./thick;
        y3_0=lcmin0.*10.^(-9).*x0C;
        y4_0=lnmin0.*10.^(-9).*x0NL;
        y5_0=rcmin0.*10.^(-9).*x0C;
        y6_0=rnmin0.*10.^(-9).*x0NR;
[T,Y]=ode23s(@thorn,tspan,[y1_0 y2_0 y3_0 y4_0 y5_0 y6_0]); % solver

% converting outputs back into input units 
OGlleaf=Y(:,1);
OGleaf=(Y(:,1)).*theta;
OGroot=(Y(:,2)).*thick;
OGleafc=Y(:,3);%10.^(9).*(Y(:,3))./x0C;
OGleafn=Y(:,4);%10.^(9).*(Y(:,4))./x0NL;
OGrootc=Y(:,5);%10.^(9).*(Y(:,5))./x0C;
OGrootn=Y(:,6);%10.^(9).*(Y(:,6))./x0NR;
N1=[N1,10.^(9).*OGleafn(end)./x0NL];
C1=[C1,10.^(9).*OGleafc(end)./x0C];

OGRGRL=(h(Y(:,3),Y(:,4)).*x0.*Yg).*(86400);
OGRGRR=(g(Y(:,6),Y(:,5)).*x0.*Yg).*(86400);
% OGnns=Y(:,4).*Y(:,1).*62.*(10.^3);
% OGnnr=Y(:,6).*Y(:,2).*62.*(10.^3);
% OGNF=((nns+nnr)./plant).*100;
OGKC=f(co2,Y(:,3),Y(:,4)).*den.*(10.^(-6)).*10.^9;
OGKN=u(n,Y(:,6),Y(:,5)).*rif.*(10.^(-6)).*10.^9;
kc1=OGKC(end);
kcc1=[kcc1,kc1];
kn1=OGKN(end);
knn1=[knn1,kn1];
% OGKN=u(n,Y(:,6),Y(:,5)).*rif.*(10.^(-6)).*10.^9;
%time=T./54; %converting seconds into months
OGtime=T./86400; %converting seconds into days
    else
        %% after defoliation
        tspann=linspace(864100,1728000,8641);
        x1_0=OGlleaf(end)./2; % halving leaf mass
        x2_0=OGroot(end)./thick;
        x3_0=OGleafc(end);%.*10.^(-9).*x0C;
        x4_0=OGleafn(end);%.*10.^(-9).*x0NL;
        x5_0=OGrootc(end);%.*10.^(-9).*x0C;
        x6_0=OGrootn(end);%.*10.^(-9).*x0NR;
[S,Q]=ode23s(@thorn,tspann,[x1_0 x2_0 x3_0 x4_0 x5_0 x6_0]); % solver
% converting outputs back into colin units 
leaf=[OGleaf ;(Q(:,1)).*theta];
root=[OGroot ;(Q(:,2)).*thick];
leafc=[OGleafc.*(10.^9)./x0C ;10.^(9).*(Q(:,3))./x0C];
leafn=[OGleafn.*(10.^9)./x0NL ;10.^(9).*(Q(:,4))./x0NL];
rootc=[OGrootc.*(10.^9)./x0C ;10.^(9).*(Q(:,5))./x0C];
rootn=[OGrootn.*(10.^9)./x0NR ;10.^(9).*(Q(:,6))./x0NR];
LAFN=10.^(9).*(Q(:,4))./x0NL;
LAFC=10.^(9).*(Q(:,3))./x0C;
% SR=Y(:,1)./Y(:,2);
% RS=Y(:,2)./Y(:,1);
 SR=leaf./root;
 RS=root./leaf;
%SR2=leaf./root;
plant=leaf+root;
carbon=leafc+rootc;
nitrogen=leafn+rootn;
ceaf=max(leafc)
coot=max(rootc);
noot=max(rootn);
neaf=max(leafn);
RGRL=[OGRGRL;(h(Q(:,3),Q(:,4)).*x0.*Yg).*(86400)];
RGRR=[OGRGRR ;(g(Q(:,6),Q(:,5)).*x0.*Yg).*(86400)];
nns=Q(:,4).*Q(:,1).*62.*(10.^3);
nnr=Q(:,6).*Q(:,2).*62.*(10.^3);
% NF=((nns+nnr)./plant).*100;
 KC=f(co2,Q(:,3),Q(:,4)).*den.*(10.^(-6)).*10.^9;
 KCC=[OGKC;KC];
 KN=u(n,Q(:,6),Q(:,5)).*rif.*(10.^(-6)).*10.^9;
 KNN=[OGKN;KN];
 RAR=S<=1469000;
 rar1=sum(RAR);
 ndkc=KC(rar1);
 kc2=KC(end);
 kcc2=[kcc2,ndkc];
 rar2=sum(RAR)
 ndkn=KN(rar2);
 knn2=[knn2,ndkn];
 C2=[C2,LAFC(rar1)];
 CUP=[C1 C2];
 N2=[N2,LAFN(rar1)];
NUP=[N1 N2];
% KN=u(n,Y(:,6),Y(:,5)).*rif.*(10.^(-6)).*10.^9;
%time=T./54; %converting seconds into months
time=[T;S]./86400; %converting seconds into days
whos time
    end
end
% creates markers at time 10 and 17 days 
cup=[kcc1 kcc2]
nup=[knn1 knn2]
tttime=[10 10 17 17];
%% Making figures
% Plot of carbon uptake rate
figure(11)
            plot(time,KCC,'color',C{i},'LineWidth',4);
            xlabel('Time (days)','FontSize',50)
            ylabel('Carbon uptake rate  (\mu molm^{-2}s^{-1})','FontSize',50)
              set(gca,'LineWidth',2,'FontSize',50)
              grid on
            hold on, drawnow        
            set(gcf, 'PaperUnits', 'inches');
          %  axis([0 20 0 0.25])
 x_width=15.65 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
% Plot of nitrogen uptake rate
 figure(1066)
            plot(time,KNN,'color',C{i},'LineWidth',4);
            xlabel('Time (days)','FontSize',50)
            ylabel('Nitrogen uptake rate  (\mu molm^{-2}s^{-1})','FontSize',50)
              set(gca,'LineWidth',2,'FontSize',50)
              grid on
            hold on, drawnow        
            set(gcf, 'PaperUnits', 'inches');
          %  axis([0 20 0 0.25])
 x_width=15.65 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 
% plot of total plant mass
figure(1)
            plot(time,leaf+root,'color',C{i},'LineWidth',4);
            xlabel('Time (days)','FontSize',50)
            ylabel('Total plant mass (g)','FontSize',50)
              set(gca,'LineWidth',2,'FontSize',50)
            hold on, drawnow
            legend('350ppm CO_2','700ppm CO_2','Feedback 2','Feedback 3','Feedback 4','Feedback 5','Feedback 6','Location','Best')
             set(gcf, 'PaperUnits', 'inches');
          %  axis([0 20 0 0.25])
 x_width=15.65 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
print('FEBDFOL1','-djpeg','-loose');

% plot of leaf nitrogen
figure(222)
plot(time,leafn,'color',C{i},'LineWidth',4);
 xlabel('Time (days)','FontSize',50)
            ylabel('Leaf nitrogen nmolmg^{-1} ','FontSize',50)
              set(gca,'LineWidth',2,'FontSize',50)
            hold on, drawnow
%axis([0 40 0 35000])
             set(gcf, 'PaperUnits', 'inches');
 x_width=15.65 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
% plot of leaf carbon
figure(212)
plot(time,leafc,'color',C{i},'LineWidth',4);
 xlabel('Time (days)','FontSize',50)
            ylabel('Leaf carbon nmolmg^{-1} ','FontSize',50)
              set(gca,'LineWidth',2,'FontSize',50)
            hold on, drawnow
%axis([0 40 0 35000])
             set(gcf, 'PaperUnits', 'inches');
 x_width=15.65 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);


% Plotting concentrations over time
figure(2)
subplot(1,2,1)
plot(time,leafc,'color',C{i},'marker','square','MarkerSize',3,'LineWidth',4), xlabel('Time (days)','FontSize',45),ylabel('Concentration nmolmg^{-1}','FontSize',45); title('A');
set(gca,'LineWidth',2,'FontSize',45)
hold on, drawnow
plot(time,leafn,'color',C{i},'LineWidth',4);
legend('Leaf C 0','Leaf N 0','Leaf C 1','Leaf N 1','Leaf C 2','Leaf N 2','Leaf C 3','Leaf N 3','Leaf C 4','Leaf N 4','Leaf C 5','Leaf N 5','Leaf C 6','Leaf N 6','Location','Best');
subplot(1,2,2)
plot(time,rootc,'color',C{i},'marker','square','MarkerSize',3,'LineWidth',4), xlabel('Time (days)','FontSize',45),ylabel('Concentration nmolmg^{-1}','FontSize',45);title('B');
set(gca,'LineWidth',2,'FontSize',45)
hold on, drawnow
plot(time,rootn,'color',C{i},'LineWidth',4);
legend('Root C 0','Root N 0','Root C 1','Root N 1','Root C 2','Root N 2','Root C 3','Root N 3','Root C 4','Root N 4','Root C 5','Root N 5','Root C 6','Root N 6','Location','Best');
hold on, drawnow
 set(gcf, 'PaperUnits', 'inches');
 x_width=31.3 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
print('FEBDFOL3','-djpeg','-loose');

% Plot of shoot:root ratio
figure(4)
%plot(time,SR,'color',C{i},'LineWidth',4);
%hold on, drawnow
plot(time,SR,'color',C{i},'LineWidth',4);
xlabel('Time (days)','FontSize',50);
ylabel('Shoot:root ratio','FontSize',50);
 set(gca,'LineWidth',2,'FontSize',50);
 hold on, drawnow 
% axis([0 40 0.3 1])
  legend('350ppm CO_2','700ppm CO_2','Location','Best');
 set(gcf, 'PaperUnits', 'inches');
 x_width=15.65 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
print('FEBDFOL4','-djpeg','-loose'); 
% plotting leaf mass and root mass over time
  figure(3) 
plot(time,leaf,'color',C{i},'marker','square','MarkerSize',3,'LineWidth',4);
xlabel('Time (days)','FontSize',50);
ylabel('Mass (g)','FontSize',50);
 set(gca,'LineWidth',2,'FontSize',50);
 hold on, drawnow 
plot(time,root,'color',C{i},'LineWidth',4);
hold on, drawnow
legend('Leaf 350ppm ','Root 350ppm ','Leaf 700ppm','Root 700ppm','Location','Best');
grid on 
 set(gcf, 'PaperUnits', 'inches');
 x_width=15.65 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
print('FEBDFOL5','-djpeg','-loose');

 % plotting growth rate change over time
figure(5) 
RGR=gradient(log(plant));
L1=length(plant);
L2=length(RGR);
sam=max(RGR);
fajitas=max(fa);
burrito=min(fa);
%mcarb=max(carbon);
%mnit=max(nitrogen)
%fa=[fa,mcarb];
%fb=[fb,mnit];
plot(time,RGRL,'color',C{i},'marker','square','MarkerSize',3,'LineWidth',4);
hold on, drawnow
plot(time,RGRR,'color',C{i},'LineWidth',4);
xlabel('Time (days)','FontSize',50)
ylabel('Relative growth rate','FontSize',50)
set(gca,'LineWidth',2,'FontSize',50)
%axis([0 40 0 0.12])
  hold on, drawnow
 set(gcf, 'PaperUnits', 'inches');
 legend('Leaf 350ppm ','Root 350ppm ','Leaf 700ppm ','Root 700ppm ','Location','Best');
 x_width=15.65 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
print('FEBDFOL6','-djpeg','-loose');
leaf(end)+root(end)
end
 % code to add the markers to the graphs
 figure(11)
        plot(tttime,cup,'ko','MarkerSize',15,'MarkerFaceColor',[0.25, 0.25, 0.25]);
                    legend('350ppm CO_2','700ppm CO_2','Location','Best')
                        print('nof2KCDFOL1','-djpeg','-loose');
 figure(1066)
        plot(tttime,nup,'ko','MarkerSize',15,'MarkerFaceColor',[0.25, 0.25, 0.25]);
                    legend('350ppm CO_2','700ppm CO_2','Location','Best')
                  print('KNDFOL1','-djpeg','-loose');

figure(222)                   
        plot(tttime,NUP,'ko','MarkerSize',15,'MarkerFaceColor',[0.25, 0.25, 0.25]);
legend('350ppm CO_2','700ppm CO_2','Location','Best')
print('nof2FEBDFOL2','-djpeg','-loose');

figure(212)
        plot(tttime,CUP,'ko','MarkerSize',15,'MarkerFaceColor',[0.25, 0.25, 0.25]);
        legend('350ppm CO_2','700ppm CO_2','Location','Best')
print('nof2FEBDFOL3','-djpeg','-loose');
%% Defining the functions for uptake rates, growth rates and model ODEs
    % Carbon uptake rate
    function ff=f(z,x,y)
        ff=(Kc.*z)./(z+kcc)-x.*a1.*aa+((((Kc.*z)./(z+kcc))./4)./(1+100000.*exp(-100.*(y-j)))).*ee;
    end
% Nitrogen uptake rate
    function uu=u(z,x,y)
        uu=(Kn.*z)./(z+knn)-x.*burger.*bb+(((((Kn.*z)./(z+knn))./4)./(1+100000.*exp(-100.*((y-j))))).*f1);
    end
%  6 Model ODEs
    function dydt=thorn(~,y)
       B=(y(1)+y(2));
        dydt=zeros(6,1);
        dydt(1)= x0.*Yg.*y(1).*h(y(3),y(4)); % eq 1: production of new shoot material (Vs)
        dydt(2)= x0.*Yg.*y(2).*g(y(6),y(5)); % eq 2:production of new root material (Vr) 
        dydt(3)= f(co2,y(3),y(4))-B.*((y(3)-y(5))./(rc.*y(1)))-(h(y(3),y(4)).*(1+y(3).*x0.*Yg))-R1; % eq 3: change in shoot carbon (Cs)
        dydt(4)= B.*((y(6)-y(4))./(rn.*y(1)))-(LL+y(4).*x0.*Yg).*h(y(3),y(4)); % eq 4: change in shoot nitrogen (Ns)
        dydt(5)= B.*((y(3)-y(5))./(rc.*y(2)))-(1+y(5).*x0.*Yg).*g(y(6),y(5))-R2  ; % eq 5: change in root carbon (Cr)
        dydt(6)= u(n,y(6),y(5))-B.*((y(6)-y(4))./(rn.*y(2)))-(LR+y(6).*x0.*Yg).*g(y(6),y(5))   ; % eq 6: change in root nitrogen (Nr)
    end
% root growth rate 
    function gg=g(x,y)
        gg=(((V1.*x)./(x+k1)).*((V2.*y)./(k2+y))).*(y.^(d1));
    end
% leaf growth rate 
 function hh=h(x,y)
        hh=(((V1.*x)./(x+k1)).*((V2.*y)./(k2+y))).*(y.^cc);
    end
 
end
