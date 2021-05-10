function co2Ntest
% Code for solving thornley model converted thornley's input values 
close ALL
fa=[];
fb=[];
alun=[];



%% input varibles
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
lcmax0=1997.6;
lcmin0=92.8;
lnmax0=1358.4;
lnmin0=0.1;
rcmax0=392.7;
rcmin0=63;
rnmax0=996.4;
rnmin0=7.54;
k1=1000; k2=1000;  % k1,k2 and v1,v2 control RGR 
V1=60;
V2=60;  
r1=(15).*10.^(-12);% leaf respiration rate kgmolg^-1s^-1
r2=(10).*10.^(-12);% root respiration rate kgmolg^-1s^-1
ld=10.^6; %leaf tissue density gm^-3
rd=0.2.*10.^6; % root tissue density gm^-3
R1=r1.*ld
R2=r2.*thick

knn=103;
kcc=200;
%% Converting to thornley units/dimensions 
%-------------------------------------------------------------------------
% definitions
x0=(1./(theta)).*((10.^3).*12) % "theta" conversion factor for plant material to volume converting input into thornley units: m^3*(kgmol)^-1
rc=0.5.*10.^(3); % carbon transport resistence
rn=1.*10.^(3); % nitrogen transport resistence
Kn=sky.*thick.*(10.^(-9).*10.^(-3)).*1; % Nitrogen uptake rate converting micromol*kg^(-1)*s^(-1) to kg*m^-3 by 10^-9 then 10^-9 from micromol to kgmol. thornley units: kgmol*m^-3s^-1
Kc=(10.^(-3).*fly./(den));%./4.82; % Carbon uptake rate=photosynthesis/(leaf thickness in metres) converted from mumol to kgmol
x0C=theta.*0.45;%Theta conversion for c concentration
x0NL=theta.*0.45.*LL;% theta conversion for n concentration
x0NR=theta.*0.45.*LR;
%x1=0.3;
ileaf=0.01;
iroot=0.01;
j=400.*x0C.*10.^(-9);
a1=10.^(-5);burger=10.^(-5);
%burger=0.45;
lam=Kc./Kn
% converting leaf concentrations into thornley units
lcmax=lcmax0.*x0C.*10.^(-9)
lcmin=lcmin0.*x0C.*10.^(-9)
lnmax=lnmax0.*x0NL.*10.^(-9)
lnmin=lnmin0.*x0NL.*10.^(-9)
% converting root concentrations into thornley units
rcmax=rcmax0.*x0C.*10.^(-9)
rcmin=rcmin0.*x0C.*10.^(-9)
rnmax=rnmax0.*x0NR.*10.^(-9)
rnmin=rnmin0.*x0NR.*10.^(-9)

tspan=[0:100:1728000].*2; %time span for 40 days 
%tspan=[0:200:4665600]; % time span for 2years
% initial conditions 
y1_0=ileaf./theta; % initial leaf mass
y2_0=iroot./thick; % initial root mass
y3_0=lcmin0.*10.^(-9).*x0C;  % initial shoot carbon concentration
y4_0=lnmin0.*10.^(-9).*x0NL; % initial shoot nitrogen concentration
y5_0=rcmin0.*10.^(-9).*x0C;  % initial root carbon concentration
y6_0=rnmin0.*10.^(-9).*x0NR; % initial root nitrogen concentration
 %% Running the model
 C={'r','b','g','k','y','c','m'}
    % switching the feedbacks on 
             aa=1;bb=1;cc=1;d1=1;ee=1;f1=1; 

n=400;
     for i=1:2
         aac=[];
        cii=[];
         if i==1
            Aco2=350; 
         else i==2
        Aco2=700;         
         end
%Aco2=400; %atmospheric co2
co2=0.7.*Aco2; % intercellular co2 
%n=400;
% Solver
%options = odeset('RelTol',1e-11,'AbsTol',1e-13);
[T,Y]=ode23tb(@thorn,tspan,[y1_0 y2_0 y3_0 y4_0 y5_0 y6_0]); % solver

%% converting outputs back into input units 
leaf=(Y(:,1)).*theta;
root=(Y(:,2)).*thick;
leafc=10.^(9).*(Y(:,3))./x0C;
leafn=10.^(9).*(Y(:,4))./x0NL;
rootc=10.^(9).*(Y(:,5))./x0C;
rootn=10.^(9).*(Y(:,6))./x0NR;
% SR=Y(:,1)./Y(:,2);
% RS=Y(:,2)./Y(:,1);
 SR=leaf./root; % shoot:root ratio 
 RS=root./leaf;
%SR2=leaf./root;
plant=leaf+root;
carbon=leafc+rootc;
nitrogen=leafn+rootn;
ceaf=max(leafc)
coot=max(rootc);
noot=max(rootn);
neaf=max(leafn);
RGRL=(h(Y(:,3),Y(:,4)).*x0.*Yg).*(86400); % leaf RGR
RGRR=(g(Y(:,6),Y(:,5)).*x0.*Yg).*(86400); % root RGR
nns=Y(:,4).*Y(:,1).*62.*(10.^3);
nnr=Y(:,6).*Y(:,2).*62.*(10.^3);
NF=((nns+nnr)./plant).*100; % percentage of N in whole plant 
KC=f(co2,Y(:,3),Y(:,4)).*den.*(10.^(-6)).*10.^9; % carbon uptake rate
KN=u(n,Y(:,6),Y(:,5)).*rif.*(10.^(-6)).*10.^9;   % nitrogen uptake rate
%time=T./54; %converting seconds into months
time=T./86400; %converting seconds into days
%Plot of total plant volume over time wth multiple simulations
scend=Y(:,3);
snend=Y(:,4);
s1=scend(end);
s2=snend(end);
% code to calculate intercellular co2 
for ci=0:1:1000
    cii=[cii,ci];
ac=f(ci,s1,s2).*den.*(10.^(-6)).*10.^9;
aac=[aac,ac];
end
%-------------------------------------------------------------------------
% plot of carbon and nitrogen uptake rates over time 
figure(999)
plot(time,KC,'color',C{i},'marker','square','MarkerSize',3,'LineWidth',4);
hold on, drawnow
plot(time,KN,'color',C{i},'LineWidth',4);
hold on,drawnow
xlabel('Time (days)','FontSize',50)
ylabel(' Nitrogen uptake rate (\mu molm^{-2}s^{-1})','FontSize',50)
  set(gcf, 'PaperUnits', 'inches');
 % axis([0 40 4 9])        
  set(gca,'LineWidth',2,'FontSize',50)
                legend('350ppm CO_2','700ppm CO_2','Location','Best')
 x_width=15.65 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 % print('NOFhinuptake1','-depsc','-loose');
 
% plot of A/ci curve
figure(55)
plot(cii,aac,'color',C{i},'LineWidth',4);
hold on, drawnow
xlabel('Intercellular CO_2 (\mu molmol^{-1})','FontSize',50)
ylabel(' Carbon uptake rate (\mu molm^{-2}s^{-1})','FontSize',50)
  set(gcf, 'PaperUnits', 'inches');
%  axis([0 40 4 9])        
  set(gca,'LineWidth',2,'FontSize',50)
                legend('350ppm CO_2','700ppm CO_2','Location','Best')
 x_width=15.65 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 print('loACi1','-depsc','-loose');
 
% plot of nitrogen uptake rate over time
figure(99)
plot(time,KN,'color',C{i},'LineWidth',4);
hold on, drawnow
xlabel('Time (days)','FontSize',50)
ylabel(' Nitrogen uptake rate (\mu molm^{-2}s^{-1})','FontSize',50)
  set(gcf, 'PaperUnits', 'inches');
 % axis([0 40 4 9])        
  set(gca,'LineWidth',2,'FontSize',50)
                legend('350ppm CO_2','700ppm CO_2','Location','Best')
 x_width=15.65 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
% print('NOFhinuptaken1','-depsc','-loose');


% plot of nitrogen percentage over time
figure(11)
plot(time,NF,'color',C{i},'LineWidth',4);
 xlabel('Time (days)','FontSize',50)
            ylabel('% of N in whole plant','FontSize',50)
              set(gca,'LineWidth',2,'FontSize',50)
            hold on, drawnow
         %     axis([0 40 0 10])
            legend('350ppm CO_2','700ppm CO_2','Feedback 2','Feedback 3','Feedback 4','Feedback 5','Feedback 6','Location','Best')
             set(gcf, 'PaperUnits', 'inches');
 x_width=15.65 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 print('loNnof2appendix9','-depsc','-loose');


 % plot of plant mass over time 
figure(1)
            plot(time,leaf+root,'color',C{i},'LineWidth',4);
            xlabel('Time (days)','FontSize',50)
            ylabel('Total plant mass (g)','FontSize',50)
              set(gca,'LineWidth',2,'FontSize',50)
            hold on, drawnow
            legend('350ppm CO_2','700ppm CO_2','Feedback 2','Feedback 3','Feedback 4','Feedback 5','Feedback 6','Location','Best')
             set(gcf, 'PaperUnits', 'inches');
             axis([0 40 0 0.28])
 x_width=15.65 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 print('NOFhinCNTES','-depsc','-loose');
 
 % plot of plant nitrogen over time
figure(222)
plot(time,leafn+rootn,'color',C{i},'LineWidth',4);
 xlabel('Time (days)','FontSize',50)
            ylabel('Total nitrogen nmolmg^{-1} ','FontSize',50)
              set(gca,'LineWidth',2,'FontSize',50)
            hold on, drawnow
%axis([0 40 0 35000])
legend('350ppm CO_2','700ppm CO_2','Feedback 2','Feedback 3','Feedback 4','Feedback 5','Feedback 6','Location','Best')
             set(gcf, 'PaperUnits', 'inches');
 x_width=15.65 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 %print('NOFhinnitrogen11','-depsc','-loose');

% plotting shoot:root ratio over time 
figure(4)
plot(time,RS,'color',C{i},'LineWidth',4);
xlabel('Time (days)','FontSize',50);
ylabel('Root:shoot ratio','FontSize',50);
 set(gca,'LineWidth',2,'FontSize',50);
 hold on, drawnow 
 %axis([0 40 0.3 1])
  legend('350ppm CO_2','700ppm CO_2','Location','Best');
 set(gcf, 'PaperUnits', 'inches');
 x_width=15.65 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
print('loNnof2appendix3','-depsc','-loose'); 


%  % plotting growth rate change over time
figure(5) 
plot(time,RGRL,'color',C{i},'marker','square','MarkerSize',3,'LineWidth',4);
hold on, drawnow
plot(time,RGRR,'color',C{i},'LineWidth',4);
xlabel('Time (days)','FontSize',50)
ylabel('Relative growth rate','FontSize',50)
set(gca,'LineWidth',2,'FontSize',50)
%axis([0 40 0 0.12])
  hold on, drawnow
 set(gcf, 'PaperUnits', 'inches');
 legend('Leaf 350ppm CO_2','Root 350ppm CO_2','Leaf 700ppm CO_2','Root 700ppm CO_2','Location','Best');
 x_width=15.65 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
% print('NOFhinnitrogen55','-depsc','-loose');

% Relative growth rates vs plant nitrogen content 
 figure(7)
plot(leafn+rootn,RGRL,'color',C{i},'marker','square','MarkerSize',3,'LineWidth',4);
hold on, drawnow
plot(leafn+rootn,RGRR,'color',C{i},'LineWidth',4);
xlabel('Nitrogen nmolmg^{-1}','FontSize',50)
ylabel('Relative growth rate','FontSize',50)
set(gca,'LineWidth',2,'FontSize',50)
  hold on, drawnow
 set(gcf, 'PaperUnits', 'inches');
 legend('Leaf 350ppm CO_2','Root 350ppm CO_2','Leaf 700ppm CO_2','Root 700ppm CO_2','Location','Best');
 x_width=15.65 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
%print('NOFhinnitrogen77','-depsc','-loose');

% Plotting leaf concentrations over time
figure(2891)
plot(time,leafc,'color',C{i},'marker','square','MarkerSize',3,'LineWidth',4), xlabel('Time (days)','FontSize',50),ylabel('Concentration nmolmg^{-1}','FontSize',50); 
hold on, drawnow
plot(time,leafn,'color',C{i},'LineWidth',4);
hold on,drawnow
legend('Leaf C 0','Leaf N 0','Leaf C 1','Leaf N 1','Leaf C 2','Leaf N 2','Leaf C 3','Leaf N 3','Leaf C 4','Leaf N 4','Leaf C 5','Leaf N 5','Leaf C 6','Leaf N 6','Location','Best');
set(gca,'LineWidth',2,'FontSize',50);
 set(gcf, 'PaperUnits', 'inches');
 x_width=15.65 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
%print('NOFhinnitrogen8','-depsc','-loose');

% root concentrations over time 
figure(333333)
plot(time,rootc,'color',C{i},'marker','square','MarkerSize',3,'LineWidth',4), xlabel('Time (days)','FontSize',50),ylabel('Concentration nmolmg^{-1}','FontSize',50);
set(gca,'LineWidth',2,'FontSize',50)
hold on, drawnow
plot(time,rootn,'color',C{i},'LineWidth',4);
legend('Root C 0','Root N 0','Root C 1','Root N 1','Root C 2','Root N 2','Root C 3','Root N 3','Root C 4','Root N 4','Root C 5','Root N 5','Root C 6','Root N 6','Location','Best');
hold on, drawnow
 set(gcf, 'PaperUnits', 'inches');
 x_width=15.65 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
 %print('NOFhinnitogen9','-depsc','-loose');
 fpm=leaf(end)+root(end)
     end
     %% Defining the functions for uptake rates, growth rates and system of equations for model
    % Carbon uptake rate
    function ff=f(z,x,y)
        ff=(Kc.*z)./(z+kcc)-x.*a1.*aa+((((Kc.*z)./(z+kcc))./4)./(1+100000.*exp(-100.*(y-j)))).*ee;
    end
% Nitrogen uptake rate 
    function uu=u(z,x,y)
        uu=(Kn.*z)./(z+knn)-x.*burger.*bb+(((((Kn.*z)./(z+knn))./4)./(1+100000.*exp(-100.*((y-j))))).*f1);
    end
%  6 model ODEs
    function dydt=thorn(~,y)
       B=(y(1)+y(2));
        dydt=zeros(6,1);
        dydt(1)=x0.*Yg.*y(1).*h(y(3),y(4)); % eq 1: production of new shoot material (Vs)
        dydt(2)=x0.*Yg.*y(2).*g(y(6),y(5)); % eq 2:production of new root material (Vr) 
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
