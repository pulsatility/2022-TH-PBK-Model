clear all
clc
tic
%--------- Time unit: second, concentration unit: pM, volume unit: L-------%


%% -------------------------- Default Parameters ---------------------------- %%
BW    = 75;
QTC   = 0.015;
QLC   = 0.25;
QRBC  = 1 - QTC - QLC;
VTC   = 0.0003;
VLC   = 0.0257;
VPC   = 0.0424;
VRBC  = 0.91 - VTC - VLC - VPC;
HCT   = 0.44;
VTBC  = 0.18;
VLBC  = 0.11;
VRBBC = 0.0236; 
VT    = VTC*BW;
VRB   = VRBC*BW;
VL    = VLC*BW;
VPtot = VPC*BW;  %Total plasma volume 

param.QC   = 15*BW^0.74/3600*(1-HCT);  %Plasma cardiac output
param.QT   = QTC*param.QC;  %Plasma flow rate to Thyroid
param.QRB  = QRBC*param.QC;  %Plasma flow rate to RB
param.QL   = QLC*param.QC;  %Plasma flow rate to Liver
param.VTB  = VT*VTBC*(1-HCT);  %Plasma volume in Thyroid
param.VRBB = VRB*VRBBC*(1-HCT);  %Plasma volume in RB
param.VRBT = VRB*(1-VRBBC);  %Tissue volume in RB
param.VLB  = VL*VLBC*(1-HCT);  %Plasma volume in Liver
param.VLT  = VL*(1-VLBC);  %Tissue volume in Liver
param.VB   = VPtot - param.VTB - param.VRBB - param.VLB;  %Plasma volume of Body Blood compartment

param.fuT4RBT = 0.1;    
param.fuT3RBT = 0.01;   
param.fuT4LT  = 0.1;
param.fuT3LT  = 0.01;

param.TBGtot = 3.515E5*VPtot;  
param.TTRtot = 5.35E6*VPtot;   
param.ALBtot = 6.45E8*VPtot;    

param.kdT4TBG = 60;
param.kdT4TTR = 5000;
param.kdT4ALB = 1.33E6;
param.kdT3TBG = 1100;
param.kdT3TTR = 3.25E5;
param.kdT3ALB = 9.75E6;

param.k2 = 0.018;  
param.k4 = 0.0832; 
param.k6 = 1.3; 
param.k8 = 0.165; 
param.k10 = 0.69; 
param.k12 = 2.2; 

param.k1 = param.k2/param.kdT4TBG;
param.k3 = param.k4/param.kdT4TTR;
param.k5 = param.k6/param.kdT4ALB;
param.k7 = param.k8/param.kdT3TBG;
param.k9 = param.k10/param.kdT3TTR;
param.k11 = param.k12/param.kdT3ALB;

param.k20 = 1.6;  %pmol/S
ThyroidT4dailyprod = param.k20/1E12*776.87*1E6*24*3600;  %ug/day. T4 MW=776.87.
param.k21 = 1;
param.k22 = 1.6/14;  %pmol/S
ThyroidT3dailyprod = param.k22/1E12*651*1E6*24*3600;  %ug/day. %T3 MW=651.
param.k23 = 0.15;

a1 = 0.25; 
a2 = 0.25;

param.k24 = 1.7084E-6/param.fuT4RBT*a2;
param.k25 = 4.243;
param.k26 = 1.44E-6/param.fuT4LT*a1; 
param.k27 = 1.398;
param.k28 = 0.00141/param.fuT4RBT;
param.k29 = 8.9E-4/param.fuT3RBT;
param.k30 = 2.786E-4/param.fuT4LT; 
param.k31 = 1.9688E-3/param.fuT3LT; 
param.k32 = 1.7084E-6/param.fuT4RBT*(1-a2);
param.k33 = 6.3349E-6/param.fuT3RBT; 
param.k34 = 1.44E-6/param.fuT4LT*(1-a1); 
param.k35 = 3.4229E-5/param.fuT3LT; 

%Parameters for X production and clerance in the Body Blood compartment(VB)
param.k37 = 0;
param.k36 = 0; 

%Parameters for X binding to TTR
param.kdXTTR = 5000;
param.k39 = 0;%0.0832;
param.k38 = param.k39/param.kdXTTR;

%Parameters for X binding to TBG
param.kdXTBG = 60;
param.k43 = 0;%0.018;
param.k42 = param.k43/param.kdXTBG;

default_param = param;


%% -------------------------- Initial Condition ---------------------------- %%
init.tfT4B   = 0;
init.tT4TBGB = 0;
init.tT4TTRB = 0;
init.tT4ALBB = 0;
init.fT4B    = 15;
init.T4TBGB  = 0;
init.T4TTRB  = 0;
init.T4ALBB  = 0;
init.tfT3B   = 0;
init.tT3TBGB = 0;
init.tT3TTRB = 0;
init.tT3ALBB = 0;
init.fT3B    = 5;
init.T3TBGB  = 0;
init.T3TTRB  = 0;
init.T3ALBB  = 0;
init.TBGB    = param.TBGtot/param.VB;
init.TTRB    = param.TTRtot/param.VB;
init.ALBB    = param.ALBtot/param.VB;

init.tfT4T   = 0;
init.tT4TBGT = 0;
init.tT4TTRT = 0;
init.tT4ALBT = 0;
init.fT4T    = 15;
init.T4TBGT  = 0;
init.T4TTRT  = 0;
init.T4ALBT  = 0;
init.tfT3T   = 0;
init.tT3TBGT = 0;
init.tT3TTRT = 0;
init.tT3ALBT = 0;
init.fT3T    = 5;
init.T3TBGT  = 0;
init.T3TTRT  = 0;
init.T3ALBT  = 0;
init.TBGT    = 0;
init.TTRT    = 0;
init.ALBT    = 0;

init.tfT4RB   = 0;
init.tT4TBGRB = 0;
init.tT4TTRRB = 0;
init.tT4ALBRB = 0;
init.fT4RB    = 15;
init.T4TBGRB  = 0;
init.T4TTRRB  = 0;
init.T4ALBRB  = 0;
init.tfT3RB   = 0;
init.tT3TBGRB = 0;
init.tT3TTRRB = 0;
init.tT3ALBRB = 0;
init.fT3RB    = 5;
init.T3TBGRB  = 0;
init.T3TTRRB  = 0;
init.T3ALBRB  = 0;
init.TBGRB    = 0;
init.TTRRB    = 0;
init.ALBRB    = 0;
init.tT4RBT    = 0;
init.tT3RBT    = 0;
init.T4RBT    = 0;
init.T3RBT    = 0;

init.tfT4L   = 0;
init.tT4TBGL = 0;
init.tT4TTRL = 0;
init.tT4ALBL = 0;
init.fT4L    = 15;
init.T4TBGL  = 0;
init.T4TTRL  = 0;
init.T4ALBL  = 0;
init.tfT3L   = 0;
init.tT3TBGL = 0;
init.tT3TTRL = 0;
init.tT3ALBL = 0;
init.fT3L    = 5;
init.T3TBGL  = 0;
init.T3TTRL  = 0;
init.T3ALBL  = 0;
init.tT4LT    = 0;
init.tT3LT    = 0;
init.T4LT    = 0;
init.T3LT    = 0;

init.XB      = 0;
init.XTTRB   = 0;
init.XT      = 0;
init.XTTRT   = 0;
init.XRB     = 0;
init.XTTRRB  = 0;
init.XL      = 0;
init.XTTRL   = 0;
init.XTBGB   = 0;
init.XTBGT   = 0;
init.XTBGRB  = 0;
init.XTBGL   = 0;


%% -------------------------- Run to Steady State---------------------------- %%
y0 = cell2mat(struct2cell(init));
tspan = [0:10000:1000*24*3600]; %Running for 1000 days
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t,y] = ode15s('TH_PBK_Nonspatial_ODE_tracer',tspan, y0, options, param);

% 


%% -------------------------- Response to Tracer Addition ---------------------- %%
y0 = y(end, :);  %Using the steady-state values from the section above as the initial values
param = default_param;

% Adding T4 tracer or T3 tracer to Body Blood compartment
y0(1) = 45; %tfT4B 
%y0(9) = 5; %tfT3B

% Varying parameter values
%%Figs. S1-S2
% param.QC = 1 * default_param.QC; %change it to 0.5x, 1x, or 2x 

%%Fig. 3A-3C
% param.k25 = 1 * default_param.k25; %change it to 0.5x, 1x, or 2x 
% param.k39 = 1 * default_param.k30; %change it to 0.5x, 1x, or 2x 

%%Fig. 3D-3F
% param.k21 = 1 * default_param.k21; %change it to 0.5x, 1x, or 2x 
% param.k28 = 1 * default_param.k28; %change it to 0.5x, 1x, or 2x 

%%Fig. S3A-S3C
% param.k27 = 1 * default_param.k27; %change it to 0.5x, 1x, or 2x 
% param.k31 = 1 * default_param.k27; %change it to 0.5x, 1x, or 2x 

%%Fig. S3D-S3F
% param.k23 = 1 * default_param.k23; %change it to 0.5x, 1x, or 2x 
% param.k29 = 1 * default_param.k29; %change it to 0.5x, 1x, or 2x 

%%Figs. 4 and S4. Comment or uncomment 3 lines below to achaieve different combinations of THBP availability
% param.TBGtot = 0;  
% param.TTRtot = 0;   
% param.ALBtot = 0; 

%%Fig. S7A-S7C
% param.k1 = 1 * default_param.k1; %change it to 0.1x, 1x, or 10x 
% param.k2 = 1 * default_param.k2; %change it to 0.1x, 1x, or 10x 
% param.k3 = 1 * default_param.k3; %change it to 0.1x, 1x, or 10x 
% param.k4 = 1 * default_param.k4; %change it to 0.1x, 1x, or 10x 
% param.k5 = 1 * default_param.k5; %change it to 0.1x, 1x, or 10x 
% param.k6 = 1 * default_param.k6; %change it to 0.1x, 1x, or 10x 

%%Fig. S7D-S7F
% param.k7 = 1 * default_param.k7; %change it to 0.1x, 1x, or 10x 
% param.k8 = 1 * default_param.k8; %change it to 0.1x, 1x, or 10x 
% param.k9 = 1 * default_param.k9; %change it to 0.1x, 1x, or 10x 
% param.k10 = 1 * default_param.k10; %change it to 0.1x, 1x, or 10x 
% param.k11 = 1 * default_param.k11; %change it to 0.1x, 1x, or 10x 
% param.k12 = 1 * default_param.k12; %change it to 0.1x, 1x, or 10x 

tspan1 = [0:100:10*24*3600];

options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t1,y1] = ode15s('TH_PBK_Nonspatial_ODE_tracer',tspan1, y0, options, param);

%Simulation results
model.tfT4B   = y1(:,1);
model.tT4TBGB  = y1(:,2);
model.tT4TTRB  = y1(:,3);
model.tT4ALBB  = y1(:,4);
model.fT4B    = y1(:,5);
model.T4TBGB  = y1(:,6);
model.T4TTRB  = y1(:,7);
model.T4ALBB  = y1(:,8);
model.tfT3B    = y1(:,9);
model.tT3TBGB  = y1(:,10);
model.tT3TTRB  = y1(:,11);
model.tT3ALBB  = y1(:,12);
model.fT3B    = y1(:,13);
model.T3TBGB  = y1(:,14);
model.T3TTRB  = y1(:,15);
model.T3ALBB  = y1(:,16);
model.TBGB    = y1(:,17);
model.TTRB    = y1(:,18);
model.ALBB    = y1(:,19);

model.tfT4T   = y1(:,20);
model.tT4TBGT  = y1(:,21);
model.tT4TTRT  = y1(:,22);
model.tT4ALBT  = y1(:,23);
model.fT4T    = y1(:,24);
model.T4TBGT  = y1(:,25);
model.T4TTRT  = y1(:,26);
model.T4ALBT  = y1(:,27);
model.tfT3T    = y1(:,28);
model.tT3TBGT  = y1(:,29);
model.tT3TTRT  = y1(:,30);
model.tT3ALBT  = y1(:,31);
model.fT3T    = y1(:,32);
model.T3TBGT  = y1(:,33);
model.T3TTRT  = y1(:,34);
model.T3ALBT  = y1(:,35);
model.TBGT    = y1(:,36);
model.TTRT    = y1(:,37);
model.ALBT    = y1(:,38);

model.tfT4RB   = y1(:,39);
model.tT4TBGRB  = y1(:,40);
model.tT4TTRRB  = y1(:,41);
model.tT4ALBRB  = y1(:,42);
model.fT4RB    = y1(:,43);
model.T4TBGRB  = y1(:,44);
model.T4TTRRB  = y1(:,45);
model.T4ALBRB  = y1(:,46);
model.tfT3RB    = y1(:,47);
model.tT3TBGRB  = y1(:,48);
model.tT3TTRRB  = y1(:,49);
model.tT3ALBRB  = y1(:,50);
model.fT3RB    = y1(:,51);
model.T3TBGRB  = y1(:,52);
model.T3TTRRB  = y1(:,53);
model.T3ALBRB  = y1(:,54);
model.TBGRB    = y1(:,55);
model.TTRRB    = y1(:,56);
model.ALBRB    = y1(:,57);
model.tT4RBT    = y1(:,58);
model.tT3RBT    = y1(:,59);
model.T4RBT    = y1(:,60);
model.T3RBT    = y1(:,61);

model.tfT4L   = y1(:,62);
model.tT4TBGL  = y1(:,63);
model.tT4TTRL  = y1(:,64);
model.tT4ALBL  = y1(:,65);
model.fT4L    = y1(:,66);
model.T4TBGL  = y1(:,67);
model.T4TTRL  = y1(:,68);
model.T4ALBL  = y1(:,69);
model.tfT3L    = y1(:,70);
model.tT3TBGL  = y1(:,71);
model.tT3TTRL  = y1(:,72);
model.tT3ALBL  = y1(:,73);
model.fT3L    = y1(:,74);
model.T3TBGL  = y1(:,75);
model.T3TTRL  = y1(:,76);
model.T3ALBL  = y1(:,77);
model.tT4LT    = y1(:,78);
model.tT3LT    = y1(:,79);
model.T4LT    = y1(:,80);
model.T3LT    = y1(:,81);

model.XB      = y1(:,82);
model.XTTRB   = y1(:,83);
model.XT      = y1(:,84);
model.XTTRT   = y1(:,85);
model.XRB     = y1(:,86);
model.XTTRRB  = y1(:,87);
model.XL      = y1(:,88);
model.XTTRL   = y1(:,89);
model.XTBGB   = y1(:,90);
model.XTBGT   = y1(:,91);
model.XTBGRB  = y1(:,92);
model.XTBGL   = y1(:,93);

model.TBGL = (param.TBGtot - (model.TBGB + model.T4TBGB + model.T3TBGB + model.tT4TBGB + model.tT3TBGB + model.XTBGB)*param.VB - (model.TBGT + model.T4TBGT + model.T3TBGT + model.tT4TBGT + model.tT3TBGT + model.XTBGT)*param.VTB - (model.TBGRB + model.T4TBGRB + model.T3TBGRB + model.tT4TBGRB + model.tT3TBGRB + model.XTBGRB)*param.VRBB - (model.T4TBGL + model.T3TBGL + model.tT4TBGL + model.tT3TBGL + model.XTBGL)*param.VLB)/param.VLB;
model.TTRL = (param.TTRtot - (model.TTRB + model.T4TTRB + model.T3TTRB + model.tT4TTRB + model.tT3TTRB + model.XTTRB)*param.VB - (model.TTRT + model.T4TTRT + model.T3TTRT + model.tT4TTRT + model.tT3TTRT + model.XTTRT)*param.VTB - (model.TTRRB + model.T4TTRRB + model.T3TTRRB + model.tT4TTRRB + model.tT3TTRRB + model.XTTRRB)*param.VRBB - (model.T4TTRL + model.T3TTRL + model.tT4TTRL + model.tT3TTRL + model.XTTRL)*param.VLB)/param.VLB;
model.ALBL = (param.ALBtot - (model.ALBB + model.T4ALBB + model.T3ALBB + model.tT4ALBB + model.tT3ALBB)*param.VB - (model.ALBT + model.T4ALBT + model.T3ALBT + model.tT4ALBT + model.tT3ALBT)*param.VTB - (model.ALBRB + model.T4ALBRB + model.T3ALBRB + model.tT4ALBRB + model.tT3ALBRB)*param.VRBB - (model.T4ALBL + model.T3ALBL + model.tT4ALBL + model.tT3ALBL)*param.VLB)/param.VLB;


% -------------------------- Plot Results ----------------------------------- %
day = t1/3600/24;
hour = t1/3600;

%%Plot T4 tracer
figure(201)
semilogy(day, model.tfT4B + model.tT4TBGB + model.tT4TTRB + model.tT4ALBB)
hold on
xlabel("Day")
ylabel("Body Blood total T4 Tracer")
pbaspect([2,1,1])
box on

figure(202)
plot(day, model.tT4LT)
hold on
xlabel("Day")
ylabel("Liver tissue T4 Tracer")
pbaspect([2,1,1])
box on

% figure(2021)
% plot(day, model.tfT4L + model.tT4TBGL + model.tT4TTRL + model.tT4ALBL)
% hold on
% xlabel('Time (Day)')
% ylabel('Liver blood total T4 Tracer')

figure(203)
plot(day, model.tT4RBT)
hold on
xlabel("Day")
ylabel("RB tissue T4 Tracer")
pbaspect([2,1,1])
box on

% figure(2031)
% plot(day, model.tfT4RB + model.tT4TBGRB + model.tT4TTRRB + model.tT4ALBRB)
% hold on
% xlabel('Time (Day)')
% ylabel('RB blood total T4 Tracer')


%%Plot T3
figure(204)
semilogy(hour, model.tfT3B + model.tT3TBGB + model.tT3TTRB + model.tT3ALBB)
hold on
xlabel("Hour")
ylabel("Body Blood total T3 Tracer")
pbaspect([2,1,1])
box on

figure(205)
plot(hour, model.tT3LT)
hold on
xlabel("Hour")
ylabel("Liver tissue T3 Tracer")
pbaspect([2,1,1])
box on

% figure(2051)
% plot(hour, model.tfT3L + model.tT3TBGL + model.tT3TTRL + model.tT3ALBL)
% hold on
% xlabel('Time (Hour)')
% ylabel('Liver blood total T3 Tracer')

figure(206)
plot(hour, model.tT3RBT)
hold on
xlabel("Hour")
ylabel("RB tissue T3 Tracer")
pbaspect([2,1,1])
box on

% figure(2061)
% plot(hour, model.tfT3RB + model.tT3TBGRB + model.tT3TTRRB + model.tT3ALBRB)
% hold on
% xlabel('Time (Hour)')
% ylabel('RB blood total T3 Tracer')


%Amount of T3 tracer in Body Blood, Liver tissue, and RB tissue (Fig. S5A and S6A)
figure(207)
semilogy(hour, (model.tfT3B + model.tT3TBGB + model.tT3TTRB + model.tT3ALBB)*param.VB)
hold on
semilogy(hour, model.tT3LT*param.VLT)
semilogy(hour, model.tT3RBT*param.VRBT)
xlabel("Hour")
ylabel("T3 Tracer amount (pmol)")
pbaspect([2,1,1])
box on

%RB:Liver ratio of T3 tracer amount (Fig. S5B and S6B)
figure(208)
plot(hour, (model.tT3RBT*param.VRBT) ./ (model.tT3LT*param.VLT))
hold on
yline(8.0874) %RB:Liver ratio of steady-state endogenous T3 amounts
xlabel("Hour")
ylabel("RB:Liver ratio of T3 tracer amount")
pbaspect([2,1,1])
box on


%%
toc
%

