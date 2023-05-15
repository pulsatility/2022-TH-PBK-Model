clear all
clc
tic
%--------- Time unit: second, concentration unit: pM, volume unit: L-------%


%% -------------------------- Default Parameters -----------------------------%%
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
param.VB   = VPtot - param.VTB - param.VRBB - param.VLB; %Plasma volume of Body Blood compartment

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

param.k20 = 1.6; %pmol/S
ThyroidT4dailyprod = param.k20/1E12*776.87*1E6*24*3600; %ug/day. T4 MW=776.87.
param.k21 = 1;
param.k22 = 1.6/14;  %pmol/S
ThyroidT3dailyprod = param.k22/1E12*651*1E6*24*3600; %ug/day. %T3 MW=651.
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

%Parameters for X production and clearance in the Body Blood compartment(VB)
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


%% -------------------------- Initial Condition -----------------------------%%
init.fT4B    = 15;
init.T4TBGB  = 0;
init.T4TTRB  = 0;
init.T4ALBB  = 0;
init.fT3B    = 5;
init.T3TBGB  = 0;
init.T3TTRB  = 0;
init.T3ALBB  = 0;
init.TBGB    = param.TBGtot/param.VB;
init.TTRB    = param.TTRtot/param.VB;
init.ALBB    = param.ALBtot/param.VB;

init.fT4T    = 15;
init.T4TBGT  = 0;
init.T4TTRT  = 0;
init.T4ALBT  = 0;
init.fT3T    = 5;
init.T3TBGT  = 0;
init.T3TTRT  = 0;
init.T3ALBT  = 0;
init.TBGT    = 0;
init.TTRT    = 0;
init.ALBT    = 0;

init.fT4RB    = 15;
init.T4TBGRB  = 0;
init.T4TTRRB  = 0;
init.T4ALBRB  = 0;
init.fT3RB    = 5;
init.T3TBGRB  = 0;
init.T3TTRRB  = 0;
init.T3ALBRB  = 0;
init.TBGRB    = 0;
init.TTRRB    = 0;
init.ALBRB    = 0;
init.T4RBT    = 0;
init.T3RBT    = 0;

init.fT4L    = 15;
init.T4TBGL  = 0;
init.T4TTRL  = 0;
init.T4ALBL  = 0;
init.fT3L    = 5;
init.T3TBGL  = 0;
init.T3TTRL  = 0;
init.T3ALBL  = 0;
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


%% -------------------------- Run to Steady State (Tables 1-5, 7) ----------------------------%%
y0 = cell2mat(struct2cell(init));
tspan0 = [0:10000:1000*24*3600]; %Running for 1000 days
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t,y] = ode15s('TH_PBK_Nonspatial_ODE',tspan0, y0, options, default_param);

%Simulation results
model.fT4B    = y(:,1);
model.T4TBGB  = y(:,2);
model.T4TTRB  = y(:,3);
model.T4ALBB  = y(:,4);

model.fT3B    = y(:,5);
model.T3TBGB  = y(:,6);
model.T3TTRB  = y(:,7);
model.T3ALBB  = y(:,8);

model.TBGB    = y(:,9);
model.TTRB    = y(:,10);
model.ALBB    = y(:,11);

model.fT4T    = y(:,12);
model.T4TBGT  = y(:,13);
model.T4TTRT  = y(:,14);
model.T4ALBT  = y(:,15);
model.fT3T    = y(:,16);
model.T3TBGT  = y(:,17);
model.T3TTRT  = y(:,18);
model.T3ALBT  = y(:,19);
model.TBGT    = y(:,20);
model.TTRT    = y(:,21);
model.ALBT    = y(:,22);

model.fT4RB    = y(:,23);  
model.T4TBGRB  = y(:,24);
model.T4TTRRB  = y(:,25);
model.T4ALBRB  = y(:,26);
model.fT3RB    = y(:,27);
model.T3TBGRB  = y(:,28);
model.T3TTRRB  = y(:,29);
model.T3ALBRB  = y(:,30);
model.TBGRB    = y(:,31);
model.TTRRB    = y(:,32);
model.ALBRB    = y(:,33);
model.T4RBT    = y(:,34);
model.T3RBT    = y(:,35);

model.fT4L    = y(:,36);
model.T4TBGL  = y(:,37);
model.T4TTRL  = y(:,38);
model.T4ALBL  = y(:,39);
model.fT3L    = y(:,40);
model.T3TBGL  = y(:,41);
model.T3TTRL  = y(:,42);
model.T3ALBL  = y(:,43);
model.T4LT    = y(:,44);
model.T3LT    = y(:,45);

model.XB      = y(:,46);
model.XTTRB   = y(:,47);
model.XT      = y(:,48);
model.XTTRT   = y(:,49);
model.XRB     = y(:,50);
model.XTTRRB  = y(:,51);
model.XL      = y(:,52);
model.XTTRL   = y(:,53);
model.XTBGB   = y(:,54);
model.XTBGT   = y(:,55);
model.XTBGRB  = y(:,56);
model.XTBGL   = y(:,57);

model.TBGL = (param.TBGtot - (model.TBGB + model.T4TBGB + model.T3TBGB + model.XTBGB)*param.VB - (model.TBGT + model.T4TBGT + model.T3TBGT + model.XTBGT)*param.VTB - (model.TBGRB + model.T4TBGRB + model.T3TBGRB + model.XTBGRB)*param.VRBB - (model.T4TBGL + model.T3TBGL + model.XTBGL)*param.VLB)/param.VLB;
model.TTRL = (param.TTRtot - (model.TTRB + model.T4TTRB + model.T3TTRB + model.XTTRB)*param.VB - (model.TTRT + model.T4TTRT + model.T3TTRT + model.XTTRT)*param.VTB - (model.TTRRB + model.T4TTRRB + model.T3TTRRB + model.XTTRRB)*param.VRBB - (model.T4TTRL + model.T3TTRL + model.XTTRL)*param.VLB)/param.VLB;
model.ALBL = (param.ALBtot - (model.ALBB + model.T4ALBB + model.T3ALBB)*param.VB - (model.ALBT + model.T4ALBT + model.T3ALBT)*param.VTB - (model.ALBRB + model.T4ALBRB + model.T3ALBRB)*param.VRBB - (model.T4ALBL + model.T3ALBL)*param.VLB)/param.VLB;


%-----------Obtaining steady state values-----------%
% Body Blood Compartment
freeT4B = model.fT4B(end);
T4TBGB = model.T4TBGB(end);
T4TTRB = model.T4TTRB(end);
T4ALBB = model.T4ALBB(end);

freeT3B = model.fT3B(end);
T3TBGB = model.T3TBGB(end);
T3TTRB = model.T3TTRB(end);
T3ALBB = model.T3ALBB(end);

TBGB = model.TBGB(end);
TTRB = model.TTRB(end);
ALBB = model.ALBB(end);

XB = model.XB(end);
XTTRB = model.XTTRB(end);
XTBGB = model.XTBGB(end);

% Thyroid blood Compartment
freeT4T = model.fT4T(end);
T4TBGT = model.T4TBGT(end);
T4TTRT = model.T4TTRT(end);
T4ALBT = model.T4ALBT(end);

freeT3T = model.fT3T(end);
T3TBGT = model.T3TBGT(end);
T3TTRT = model.T3TTRT(end);
T3ALBT = model.T3ALBT(end);

TBGT = model.TBGT(end);
TTRT = model.TTRT(end);
ALBT = model.ALBT(end);

XT = model.XT(end);
XTTRT = model.XTTRT(end);
XTBGT = model.XTBGT(end);

% RB blood Compartment
freeT4RB = model.fT4RB(end);
T4TBGRB = model.T4TBGRB(end);
T4TTRRB = model.T4TTRRB(end);
T4ALBRB = model.T4ALBRB(end);

freeT3RB = model.fT3RB(end);
T3TBGRB  = model.T3TBGRB(end);
T3TTRRB  =  model.T3TTRRB(end);
T3ALBRB  = model.T3ALBRB(end);

TBGRB = model.TBGRB(end);
TTRRB = model.TTRRB(end);
ALBRB = model.ALBRB(end);

XRB = model.XRB(end);
XTTRRB = model.XTTRRB(end);
XTBGRB = model.XTBGRB(end);

% RB tissue Compartment
T4RBT = model.T4RBT(end);
T3RBT = model.T3RBT(end);

% Liver blood Compartment
freeT4L = model.fT4L(end);
T4TBGL = model.T4TBGL(end);
T4TTRL = model.T4TTRL(end);
T4ALBL = model.T4ALBL(end);

freeT3L = model.fT3L(end);
T3TBGL = model.T3TBGL(end);
T3TTRL = model.T3TTRL(end);
T3ALBL = model.T3ALBL(end);

TBGL = model.TBGL(end);
TTRL = model.TTRL(end);
ALBL = model.ALBL(end);

XL = model.XL(end);
XTTRL = model.XTTRL(end);
XTBGL = model.XTBGL(end);

% Liver tissue Compartment
T4LT = model.T4LT(end);
T3LT = model.T3LT(end);

% Calculated metrics for T4 and T3 in Body Blood
boundT4B = T4TBGB + T4TTRB + T4ALBB;
totalT4B = boundT4B + freeT4B;
freeT4Bpercentage = freeT4B/(totalT4B)*100;
T4TBGBpercentage = T4TBGB/totalT4B*100;
T4TTRBpercentage = T4TTRB/totalT4B*100;
T4ALBBpercentage = T4ALBB/totalT4B*100;

TBGT4Bsatpercentage = T4TBGB/(TBGB+T4TBGB+T3TBGB+XTBGB)*100;
TTRT4Bsatpercentage = T4TTRB/(TTRB+T4TTRB+T3TTRB+XTTRB)*100;
ALBT4Bsatpercentage = T4ALBB/(ALBB+T4ALBB+T3ALBB)*100;

boundT3B = T3TBGB + T3TTRB + T3ALBB;
totalT3B = boundT3B + freeT3B;
freeT3Bpercentage = freeT3B/(totalT3B)*100;
T3TBGBpercentage = T3TBGB/totalT3B*100;
T3TTRBpercentage = T3TTRB/totalT3B*100;
T3ALBBpercentage = T3ALBB/totalT3B*100;

TBGT3Bsatpercentage = T3TBGB/(TBGB+T4TBGB+T3TBGB+XTBGB)*100;
TTRT3Bsatpercentage = T3TTRB/(TTRB+T4TTRB+T3TTRB+XTTRB)*100;
ALBT3Bsatpercentage = T3ALBB/(ALBB+T4ALBB+T3ALBB)*100;

TBGXBsatpercentage = XTBGB/(TBGB+T4TBGB+T3TBGB+XTBGB)*100;
TTRXBsatpercentage = XTTRB/(TTRB+T4TTRB+T3TTRB+XTTRB)*100;

TBGBsatpercentage = TBGT4Bsatpercentage + TBGT3Bsatpercentage + TBGXBsatpercentage;
TTRBsatpercentage = TTRT4Bsatpercentage + TTRT3Bsatpercentage + TTRXBsatpercentage;
ALBBsatpercentage = ALBT4Bsatpercentage + ALBT3Bsatpercentage;

% Calculated metrics for T4 and T3 in Thyroid blood
boundT4T = T4TBGT + T4TTRT + T4ALBT;
totalT4T = boundT4T + freeT4T;
freeT4Tpercentage = freeT4T/(totalT4T)*100;
T4TBGTpercentage = T4TBGT/totalT4T*100;
T4TTRTpercentage = T4TTRT/totalT4T*100;
T4ALBTpercentage = T4ALBT/totalT4T*100;

TBGT4Tsatpercentage = T4TBGT/(TBGT+T4TBGT+T3TBGT+XTBGT)*100;
TTRT4Tsatpercentage = T4TTRT/(TTRT+T4TTRT+T3TTRT+XTTRT)*100;
ALBT4Tsatpercentage = T4ALBT/(ALBT+T4ALBT+T3ALBT)*100;

boundT3T = T3TBGT + T3TTRT + T3ALBT;
totalT3T = boundT3T + freeT3T;
freeT3Tpercentage = freeT3T/(totalT3T)*100;
T3TBGTpercentage = T3TBGT/totalT3T*100;
T3TTRTpercentage = T3TTRT/totalT3T*100;
T3ALBTpercentage = T3ALBT/totalT3T*100;

TBGT3Tsatpercentage = T3TBGT/(TBGT+T4TBGT+T3TBGT+XTBGT)*100;
TTRT3Tsatpercentage = T3TTRT/(TTRT+T4TTRT+T3TTRT+XTTRT)*100;
ALBT3Tsatpercentage = T3ALBT/(ALBT+T4ALBT+T3ALBT)*100;

TBGXTsatpercentage = XTBGT/(TBGT+T4TBGT+T3TBGT+XTBGT)*100;
TTRXTsatpercentage = XTTRT/(TTRT+T4TTRT+T3TTRT+XTTRT)*100;

TBGTsatpercentage = TBGT4Tsatpercentage + TBGT3Tsatpercentage + TBGXTsatpercentage;
TTRTsatpercentage = TTRT4Tsatpercentage + TTRT3Tsatpercentage + TTRXTsatpercentage;
ALBTsatpercentage = ALBT4Tsatpercentage + ALBT3Tsatpercentage;

% Calculated metrics for T4 and T3 in RB blood
boundT4RB = T4TBGRB + T4TTRRB + T4ALBRB;
totalT4RB = boundT4RB + freeT4RB;
freeT4RBpercentage = freeT4RB/(totalT4RB)*100;
T4TBGRBpercentage = T4TBGRB/totalT4RB*100;
T4TTRRBpercentage = T4TTRRB/totalT4RB*100;
T4ALBRBpercentage = T4ALBRB/totalT4RB*100;

TBGT4RBsatpercentage = T4TBGRB/(TBGRB+T4TBGRB+T3TBGRB+XTBGRB)*100;
TTRT4RBsatpercentage = T4TTRRB/(TTRRB+T4TTRRB+T3TTRRB+XTTRRB)*100;
ALBT4RBsatpercentage = T4ALBRB/(ALBRB+T4ALBRB+T3ALBRB)*100;

boundT3RB = T3TBGRB + T3TTRRB + T3ALBRB;
totalT3RB = boundT3RB + freeT3RB;
freeT3RBpercentage = freeT3RB/(totalT3RB)*100;
T3TBGRBpercentage = T3TBGRB/totalT3RB*100;
T3TTRRBpercentage = T3TTRRB/totalT3RB*100;
T3ALBRBpercentage = T3ALBRB/totalT3RB*100;

TBGT3RBsatpercentage = T3TBGRB/(TBGRB+T4TBGRB+T3TBGRB+XTBGRB)*100;
TTRT3RBsatpercentage = T3TTRRB/(TTRRB+T4TTRRB+T3TTRRB+XTTRRB)*100;
ALBT3RBsatpercentage = T3ALBRB/(ALBRB+T4ALBRB+T3ALBRB)*100;

TBGXRBsatpercentage = XTBGRB/(TBGRB+T4TBGRB+T3TBGRB+XTBGRB)*100;
TTRXRBsatpercentage = XTTRRB/(TTRRB+T4TTRRB+T3TTRRB+XTTRRB)*100;

TBGRBsatpercentage = TBGT4RBsatpercentage + TBGT3RBsatpercentage + TBGXRBsatpercentage;
TTRRBsatpercentage = TTRT4RBsatpercentage + TTRT3RBsatpercentage + TTRXRBsatpercentage;
ALBRBsatpercentage = ALBT4RBsatpercentage + ALBT3RBsatpercentage;

T4RBTissueToPlasma_Partition = T4RBT/totalT4RB;
T3RBTissueToPlasma_Partition = T3RBT/totalT3RB;

% Calculated metrics for T4 and T3 in Liver blood;
boundT4L = T4TBGL + T4TTRL + T4ALBL;
totalT4L = boundT4L + freeT4L;
freeT4Lpercentage = freeT4L/(totalT4L)*100;
T4TBGLpercentage = T4TBGL/totalT4L*100;
T4TTRLpercentage = T4TTRL/totalT4L*100;
T4ALBLpercentage = T4ALBL/totalT4L*100;

TBGT4Lsatpercentage = T4TBGL/(TBGL+T4TBGL+T3TBGL+XTBGL)*100;
TTRT4Lsatpercentage = T4TTRL/(TTRL+T4TTRL+T3TTRL+XTTRL)*100;
ALBT4Lsatpercentage = T4ALBL/(ALBL+T4ALBL+T3ALBL)*100;

boundT3L = T3TBGL + T3TTRL + T3ALBL;
totalT3L = boundT3L + freeT3L;
freeT3Lpercentage = freeT3L/(totalT3L)*100;
T3TBGLpercentage = T3TBGL/totalT3L*100;
T3TTRLpercentage = T3TTRL/totalT3L*100;
T3ALBLpercentage = T3ALBL/totalT3L*100;

TBGT3Lsatpercentage = T3TBGL/(TBGL+T4TBGL+T3TBGL+XTBGL)*100;
TTRT3Lsatpercentage = T3TTRL/(TTRL+T4TTRL+T3TTRL+XTTRL)*100;
ALBT3Lsatpercentage = T3ALBL/(ALBL+T4ALBL+T3ALBL)*100;

TBGXLsatpercentage = XTBGL/(TBGL+T4TBGL+T3TBGL+XTBGL)*100;
TTRXLsatpercentage = XTTRL/(TTRL+T4TTRL+T3TTRL+XTTRL)*100;

TBGLsatpercentage = TBGT4Lsatpercentage + TBGT3Lsatpercentage + TBGXLsatpercentage;
TTRLsatpercentage = TTRT4Lsatpercentage + TTRT3Lsatpercentage + TTRXLsatpercentage;
ALBLsatpercentage = ALBT4Lsatpercentage + ALBT3Lsatpercentage;

T4LiverTissueToPlasma_Partition = T4LT/totalT4L;
T3LiverTissueToPlasma_Partition = T3LT/totalT3L;

% Total T4 and T3 amounts in All Blood compartments
totalT4_in_all_blood = totalT4B*param.VB + totalT4T*param.VTB + totalT4RB*param.VRBB + totalT4L*param.VLB;
totalT3_in_all_blood = totalT3B*param.VB + totalT3T*param.VTB + totalT3RB*param.VRBB + totalT3L*param.VLB;

% Total T4 and T3 amounts in RB tissue
totalT4RBT = T4RBT*param.VRBT;
totalT3RBT = T3RBT*param.VRBT;

% Total T4 and T3 amounts in Liver tissue
totalT4LT = T4LT*param.VLT;
totalT3LT = T3LT*param.VLT;

% Total T4 and T3 amounts in all extrathyroidal tissues
totalT4_in_all_tissue = totalT4RBT + totalT4LT;
totalT3_in_all_tissue = totalT3RBT + totalT3LT;

% Total T4 and T3 amounts in Entire Body
totalT4_in_entire_body = totalT4_in_all_blood + totalT4_in_all_tissue;
totalT3_in_entire_body = totalT3_in_all_blood + totalT3_in_all_tissue;

% Percentage of Total T4 and total T3 amounts in All Blood Compartments
Percent_TotalT4_in_all_blood = totalT4_in_all_blood/totalT4_in_entire_body*100;
Percent_TotalT3_in_all_blood = totalT3_in_all_blood/totalT3_in_entire_body*100;

% Percentage of Total T4 and total T3 amounts in Liver tissue
Percent_TotalT4_in_liver_tissue = totalT4LT/totalT4_in_entire_body*100;
Percent_TotalT3_in_liver_tissue = totalT3LT/totalT3_in_entire_body*100;

% Percentage of Total T4 and total T3 amounts in RB tissue
Percent_TotalT4_in_RB_tissue = totalT4RBT/totalT4_in_entire_body*100;
Percent_TotalT3_in_RB_tissue = totalT3RBT/totalT3_in_entire_body*100;

% Percentage of Total T4 and total T3 amounts in all extrathyroidal tissues
Percent_TotalT4_in_all_tissue = totalT4_in_all_tissue/totalT4_in_entire_body*100;
Percent_TotalT3_in_all_tissue = totalT3_in_all_tissue/totalT3_in_entire_body*100;

% T3 daily production rate in extrathyroidal tissues
TissueT3dailyprod = (param.k24*T4RBT*param.fuT4RBT*param.VRBT + param.k26*T4LT*param.fuT4LT*param.VLT)/1E12*651*1E6*24*3600; %ug/day. %T3 MW=651.

% T4 concentration flux (pM/S) in the Body Blood Compartment
Flux_freeT4_veins = (freeT4T*param.QT + freeT4L*param.QL + freeT4RB*param.QRB) / param.VB;
Flux_freeT4_artery = freeT4B*param.QC/param.VB;

Flux_T4TBG_veins = (T4TBGT*param.QT + T4TBGL*param.QL + T4TBGRB*param.QRB) / param.VB;
Flux_T4TBG_artery = T4TBGB*param.QC/param.VB;
Flux_T4TBG_association_blood = freeT4B*TBGB*param.k1;
Flux_T4TBG_dissociation_blood = T4TBGB*param.k2;
Expected_equilibrium_T4TBG_artery = freeT4B*TBGB / param.kdT4TBG;
Percentage_difference_to_equilibrium_T4TBG_blood = (T4TBGB - Expected_equilibrium_T4TBG_artery) / Expected_equilibrium_T4TBG_artery *100;

Flux_T4TTR_veins = (T4TTRT*param.QT + T4TTRL*param.QL + T4TTRRB*param.QRB) / param.VB;
Flux_T4TTR_artery = T4TTRB*param.QC/param.VB;
Flux_T4TTR_association_blood = freeT4B*TTRB*param.k3;
Flux_T4TTR_dissociation_blood = T4TTRB*param.k4;
Expected_equilibrium_T4TTR_artery = freeT4B*TTRB / param.kdT4TTR;
Percentage_difference_to_equilibrium_T4TTR_blood = (T4TTRB - Expected_equilibrium_T4TTR_artery) / Expected_equilibrium_T4TTR_artery *100;

Flux_T4ALB_veins = (T4ALBT*param.QT + T4ALBL*param.QL + T4ALBRB*param.QRB) / param.VB;
Flux_T4ALB_artery = T4ALBB*param.QC/param.VB;
Flux_T4ALB_association_blood = freeT4B*ALBB*param.k5;
Flux_T4ALB_dissociation_blood = T4ALBB*param.k6;
Expected_equilibrium_T4ALB_artery = freeT4B*ALBB / param.kdT4ALB;
Percentage_difference_to_equilibrium_T4ALB_blood = (T4ALBB - Expected_equilibrium_T4ALB_artery) / Expected_equilibrium_T4ALB_artery *100;

Flux_T4THBP_artery = Flux_T4TBG_artery + Flux_T4TTR_artery + Flux_T4ALB_artery;
Flux_T4THBP_veins = Flux_T4TBG_veins + Flux_T4TTR_veins + Flux_T4ALB_veins;
Flux_T4THBP_association_blood = Flux_T4TBG_association_blood + Flux_T4TTR_association_blood + Flux_T4ALB_association_blood;
Flux_T4THBP_dissociation_blood = Flux_T4TBG_dissociation_blood + Flux_T4TTR_dissociation_blood + Flux_T4ALB_dissociation_blood;

% T3 concentration flux (pM/S) in the Body Blood Compartment
Flux_freeT3_veins = (freeT3T*param.QT + freeT3L*param.QL + freeT3RB*param.QRB) / param.VB;
Flux_freeT3_artery = freeT3B*param.QC/param.VB;

Flux_T3TBG_veins = (T3TBGT*param.QT + T3TBGL*param.QL + T3TBGRB*param.QRB) / param.VB;
Flux_T3TBG_artery = T3TBGB*param.QC/param.VB;
Flux_T3TBG_association_blood = freeT3B*TBGB*param.k7;
Flux_T3TBG_dissociation_blood = T3TBGB*param.k8;
Expected_equilibrium_T3TBG_artery = freeT3B*TBGB / param.kdT3TBG;
Percentage_difference_to_equilibrium_T3TBG_blood = (T3TBGB - Expected_equilibrium_T3TBG_artery) / Expected_equilibrium_T3TBG_artery *100;

Flux_T3TTR_veins = (T3TTRT*param.QT + T3TTRL*param.QL + T3TTRRB*param.QRB) / param.VB;
Flux_T3TTR_artery = T3TTRB*param.QC/param.VB;
Flux_T3TTR_association_blood = freeT3B*TTRB*param.k9;
Flux_T3TTR_dissociation_blood = T3TTRB*param.k10;
Expected_equilibrium_T3TTR_artery = freeT3B*TTRB / param.kdT3TTR;
Percentage_difference_to_equilibrium_T3TTR_blood = (T3TTRB - Expected_equilibrium_T3TTR_artery) / Expected_equilibrium_T3TTR_artery *100;

Flux_T3ALB_veins = (T3ALBT*param.QT + T3ALBL*param.QL + T3ALBRB*param.QRB) / param.VB;
Flux_T3ALB_artery = T3ALBB*param.QC/param.VB;
Flux_T3ALB_association_blood = freeT3B*ALBB*param.k11;
Flux_T3ALB_dissociation_blood = T3ALBB*param.k12;
Expected_equilibrium_T3ALB_artery = freeT3B*ALBB / param.kdT3ALB;
Percentage_difference_to_equilibrium_T3ALB_blood = (T3ALBB - Expected_equilibrium_T3ALB_artery) / Expected_equilibrium_T3ALB_artery *100;


% T4 concentration flux (pM/S) in the Liver blood Compartment
Flux_freeT4_liver_artery = freeT4B*param.QL/param.VLB;
Flux_freeT4_liver_vein = freeT4L*param.QL / param.VLB;

Flux_T4TBG_liver_artery = T4TBGB * param.QL / param.VLB;
Flux_T4TBG_liver_vein = T4TBGL * param.QL / param.VLB;
Flux_T4TBG_association_liver = freeT4L*TBGL*param.k1;
Flux_T4TBG_dissociation_liver = T4TBGL*param.k2;
Expected_equilibrium_T4TBG_liver = freeT4L*TBGL / param.kdT4TBG;
Percentage_difference_to_equilibrium_T4TBG_liver = (T4TBGL - Expected_equilibrium_T4TBG_liver) / Expected_equilibrium_T4TBG_liver *100;

Flux_T4TTR_liver_artery = T4TTRB * param.QL / param.VLB;
Flux_T4TTR_liver_vein = T4TTRL * param.QL / param.VLB;
Flux_T4TTR_association_liver = freeT4L*TTRL*param.k3;
Flux_T4TTR_dissociation_liver = T4TTRL*param.k4;
Expected_equilibrium_T4TTR_liver = freeT4L*TTRL / param.kdT4TTR;
Percentage_difference_to_equilibrium_T4TTR_liver = (T4TTRL - Expected_equilibrium_T4TTR_liver) / Expected_equilibrium_T4TTR_liver *100;

Flux_T4ALB_liver_artery = T4ALBB * param.QL / param.VLB;
Flux_T4ALB_liver_vein = T4ALBL * param.QL / param.VLB;
Flux_T4ALB_association_liver = freeT4L*ALBL*param.k5;
Flux_T4ALB_dissociation_liver = T4ALBL*param.k6;
Expected_equilibrium_T4ALB_liver = freeT4L*ALBL / param.kdT4ALB;
Percentage_difference_to_equilibrium_T4ALB_liver = (T4ALBL - Expected_equilibrium_T4ALB_liver) / Expected_equilibrium_T4ALB_liver *100;

Flux_T4THBP_liver_artery = Flux_T4TBG_liver_artery + Flux_T4TTR_liver_artery + Flux_T4ALB_liver_artery;
Flux_T4THBP_liver_vein = Flux_T4TBG_liver_vein + Flux_T4TTR_liver_vein + Flux_T4ALB_liver_vein;
Flux_T4THBP_association_liver = Flux_T4TBG_association_liver + Flux_T4TTR_association_liver + Flux_T4ALB_association_liver;
Flux_T4THBP_dissociation_liver = Flux_T4TBG_dissociation_liver + Flux_T4TTR_dissociation_liver + Flux_T4ALB_dissociation_liver;

Flux_LiverT4_clearance = (param.k26+param.k34)*T4LT*param.fuT4LT*param.VLT/param.VLB;


% T3 concentration flux (pM/S) in the Liver blood Compartment
Flux_freeT3_liver_artery = freeT3B*param.QL/param.VLB;
Flux_freeT3_liver_vein = freeT3L*param.QL / param.VLB;

Flux_T3TBG_liver_artery = T3TBGB * param.QL / param.VLB;
Flux_T3TBG_liver_vein = T3TBGL * param.QL / param.VLB;
Flux_T3TBG_association_liver = freeT3L*TBGL*param.k7;
Flux_T3TBG_dissociation_liver = T3TBGL*param.k8;
Expected_equilibrium_T3TBG_liver = freeT3L*TBGL / param.kdT3TBG;
Percentage_difference_to_equilibrium_T3TBG_liver = (T3TBGL - Expected_equilibrium_T3TBG_liver) / Expected_equilibrium_T3TBG_liver *100;

Flux_T3TTR_liver_artery = T3TTRB * param.QL / param.VLB;
Flux_T3TTR_liver_vein = T3TTRL * param.QL / param.VLB;
Flux_T3TTR_association_liver = freeT3L*TTRL*param.k9;
Flux_T3TTR_dissociation_liver = T3TTRL*param.k10;
Expected_equilibrium_T3TTR_liver = freeT3L*TTRL / param.kdT3TTR;
Percentage_difference_to_equilibrium_T3TTR_liver = (T3TTRL - Expected_equilibrium_T3TTR_liver) / Expected_equilibrium_T3TTR_liver *100;

Flux_T3ALB_liver_artery = T3ALBB * param.QL / param.VLB;
Flux_T3ALB_liver_vein = T3ALBL * param.QL / param.VLB;
Flux_T3ALB_association_liver = freeT3L*ALBL*param.k11;
Flux_T3ALB_dissociation_liver = T3ALBL*param.k12;
Expected_equilibrium_T3ALB_liver = freeT3L*ALBL / param.kdT3ALB;
Percentage_difference_to_equilibrium_T3ALB_liver = (T3ALBL - Expected_equilibrium_T3ALB_liver) / Expected_equilibrium_T3ALB_liver *100;

Flux_T3_clearance_liver =param.k35*T3LT*param.fuT3LT*param.VLT/param.VLB;


% T4 concentration flux (pM/S) in the RB blood Compartment
Flux_freeT4_RB_artery = freeT4B*param.QRB/param.VRBB;
Flux_freeT4_RB_vein = freeT4RB*param.QRB / param.VRBB;

Flux_T4TBG_RB_artery = T4TBGB * param.QRB / param.VRBB;
Flux_T4TBG_RB_vein = T4TBGRB * param.QRB / param.VRBB;
Flux_T4TBG_association_RB = freeT4RB*TBGRB*param.k1;
Flux_T4TBG_dissociation_RB = T4TBGRB*param.k2;
Expected_equilibrium_T4TBG_RB = freeT4RB*TBGRB / param.kdT4TBG;
Percentage_difference_to_equilibrium_T4TBG_RB = (T4TBGRB - Expected_equilibrium_T4TBG_RB) / Expected_equilibrium_T4TBG_RB *100;

Flux_T4TTR_RB_artery = T4TTRB * param.QRB / param.VRBB;
Flux_T4TTR_RB_vein = T4TTRRB * param.QRB / param.VRBB;
Flux_T4TTR_association_RB = freeT4RB*TTRRB*param.k3;
Flux_T4TTR_dissociation_RB = T4TTRRB*param.k4;
Expected_equilibrium_T4TTR_RB = freeT4RB*TTRRB / param.kdT4TTR;
Percentage_difference_to_equilibrium_T4TTR_RB = (T4TTRRB - Expected_equilibrium_T4TTR_RB) / Expected_equilibrium_T4TTR_RB *100;

Flux_T4ALB_RB_artery = T4ALBB * param.QRB / param.VRBB;
Flux_T4ALB_RB_vein = T4ALBRB * param.QRB / param.VRBB;
Flux_T4ALB_association_RB = freeT4RB*ALBRB*param.k5;
Flux_T4ALB_dissociation_RB = T4ALBRB*param.k6;
Expected_equilibrium_T4ALB_RB = freeT4RB*ALBRB / param.kdT4ALB;
Percentage_difference_to_equilibrium_T4ALB_RB = (T4ALBRB - Expected_equilibrium_T4ALB_RB) / Expected_equilibrium_T4ALB_RB *100;

Flux_T4THBP_RB_artery = Flux_T4TBG_RB_artery + Flux_T4TTR_RB_artery + Flux_T4ALB_RB_artery;
Flux_T4THBP_RB_vein = Flux_T4TBG_RB_vein + Flux_T4TTR_RB_vein + Flux_T4ALB_RB_vein;
Flux_T4THBP_association_RB = Flux_T4TBG_association_RB + Flux_T4TTR_association_RB + Flux_T4ALB_association_RB;
Flux_T4THBP_dissociation_RB = Flux_T4TBG_dissociation_RB + Flux_T4TTR_dissociation_RB + Flux_T4ALB_dissociation_RB;

Flux_T4_clearance_RB = (param.k24+param.k32)*T4RBT*param.fuT4RBT*param.VRBT/param.VRBB;


% T3 concentration flux (pM/S) in the RB blood Compartment
Flux_freeT3_RB_artery = freeT3B*param.QRB/param.VRBB;
Flux_freeT3_RB_vein = freeT3RB*param.QRB / param.VRBB;

Flux_T3TBG_RB_artery = T3TBGB * param.QRB / param.VRBB;
Flux_T3TBG_RB_vein = T3TBGRB * param.QRB / param.VRBB;
Flux_T3TBG_association_RB = freeT3RB*TBGRB*param.k7;
Flux_T3TBG_dissociation_RB = T3TBGRB*param.k8;
Expected_equilibrium_T3TBG_RB = freeT3RB*TBGRB / param.kdT3TBG;
Percentage_difference_to_equilibrium_T3TBG_RB = (T3TBGRB - Expected_equilibrium_T3TBG_RB) / Expected_equilibrium_T3TBG_RB *100;

Flux_T3TTR_RB_artery = T3TTRB * param.QRB / param.VRBB;
Flux_T3TTR_RB_vein = T3TTRRB * param.QRB / param.VRBB;
Flux_T3TTR_association_RB = freeT3RB*TTRRB*param.k9;
Flux_T3TTR_dissociation_RB = T3TTRRB*param.k10;
Expected_equilibrium_T3TTR_RB = freeT3RB*TTRRB / param.kdT3TTR;
Percentage_difference_to_equilibrium_T3TTR_RB = (T3TTRRB - Expected_equilibrium_T3TTR_RB) / Expected_equilibrium_T3TTR_RB *100;

Flux_T3ALB_RB_artery = T3ALBB * param.QRB / param.VRBB;
Flux_T3ALB_RB_vein = T3ALBRB * param.QRB / param.VRBB;
Flux_T3ALB_association_RB = freeT3RB*ALBRB*param.k11;
Flux_T3ALB_dissociation_RB = T3ALBRB*param.k12;
Expected_equilibrium_T3ALB_RB = freeT3RB*ALBRB / param.kdT3ALB;
Percentage_difference_to_equilibrium_T3ALB_RB = (T3ALBRB - Expected_equilibrium_T3ALB_RB) / Expected_equilibrium_T3ALB_RB *100;

Flux_T3_clearance_RB = param.k33*T3RBT*param.fuT3RBT*param.VRBT/param.VRBB;


% Blood T4 tissue influx rate (pmol/S)
LiverT4_influx_rate = param.k25*freeT4L;
RBT4_influx_rate = param.k21*freeT4RB;
BodyT4_influx_rate = LiverT4_influx_rate + RBT4_influx_rate;
Percent_LiverT4_influx_over_liver_blood_supply = LiverT4_influx_rate/(totalT4B*param.QL);

% Blood T3 tissue influx rate (pmol/S)
LiverT3_influx_rate = param.k27*freeT3L;
RBT3_influx_rate = param.k23*freeT3RB;
BodyT3_influx_rate = LiverT3_influx_rate + RBT3_influx_rate;
Percent_LiverT3_influx_over_liver_blood_supply = LiverT3_influx_rate/(totalT3B*param.QL);

% Liver T4 arterial blood supply rate (pmol/S)
freeT4_liver_blood_supply_rate = freeT4B*param.QL;
totalT4_liver_blood_supply_rate = totalT4B*param.QL;
Percent_influx_freeT4_liver_blood_supply_rate = freeT4_liver_blood_supply_rate/LiverT4_influx_rate * 100;
Percent_influx_totalT4_liver_blood_supply_rate = totalT4_liver_blood_supply_rate/LiverT4_influx_rate  * 100;

% Liver T4 venous blood exit rate (pmol/S)
freeT4_liver_blood_exit_rate = freeT4L*param.QL;
totalT4_liver_blood_exit_rate = totalT4L*param.QL;
Percent_influx_freeT4_liver_blood_exit_rate = freeT4_liver_blood_exit_rate/LiverT4_influx_rate * 100;
Percent_influx_totalT4_liver_blood_exit_rate = totalT4_liver_blood_exit_rate/LiverT4_influx_rate  * 100;

% Liver T3 arterial blood supply rate (pmol/S)
freeT3_liver_blood_supply_rate = freeT3B*param.QL;
totalT3_liver_blood_supply_rate = totalT3B*param.QL;
Percent_influx_freeT3_liver_blood_supply_rate = freeT3_liver_blood_supply_rate/LiverT3_influx_rate * 100;
Percent_influx_totalT3_liver_blood_supply_rate = totalT3_liver_blood_supply_rate/LiverT3_influx_rate * 100;

% Liver T3 venous blood exit rate (pmol/S)
freeT3_liver_blood_exit_rate = freeT3L*param.QL;
totalT3_liver_blood_exit_rate = totalT3L*param.QL;
Percent_influx_freeT3_liver_blood_exit_rate = freeT3_liver_blood_exit_rate/LiverT3_influx_rate * 100;
Percent_influx_totalT3_liver_blood_exit_rate = totalT3_liver_blood_exit_rate/LiverT3_influx_rate  * 100;


% RB T4 arterial blood supply rate (pmol/S)
freeT4_RB_blood_supply_rate = freeT4B*param.QRB;
totalT4_RB_blood_supply_rate = totalT4B*param.QRB;
Percent_influx_freeT4_RB_blood_supply_rate = freeT4_RB_blood_supply_rate/RBT4_influx_rate * 100;
Percent_influx_totalT4_RB_blood_supply_rate = totalT4_RB_blood_supply_rate/RBT4_influx_rate  * 100;

% RB T4 venous blood exit rate (pmol/S)
freeT4_RB_blood_exit_rate = freeT4RB*param.QRB;
totalT4_RB_blood_exit_rate = totalT4RB*param.QRB;
Percent_influx_freeT4_RB_blood_exit_rate = freeT4_RB_blood_exit_rate/RBT4_influx_rate * 100;
Percent_influx_totalT4_RB_blood_exit_rate = totalT4_RB_blood_exit_rate/RBT4_influx_rate  * 100;

% RB T3 venous blood supply rate (pmol/S)
freeT3_RB_blood_supply_rate = freeT3B*param.QRB;
totalT3_RB_blood_supply_rate = totalT3B*param.QRB;
Percent_influx_freeT3_RB_blood_supply_rate = freeT3_RB_blood_supply_rate/RBT3_influx_rate * 100;
Percent_influx_totalT3_RB_blood_supply_rate = totalT3_RB_blood_supply_rate/RBT3_influx_rate * 100;

% RB T3 venous blood exit rate (pmol/S)
freeT3_RB_blood_exit_rate = freeT3RB*param.QRB;
totalT3_RB_blood_exit_rate = totalT3RB*param.QRB;
Percent_influx_freeT3_RB_blood_exit_rate = freeT3_RB_blood_exit_rate/RBT3_influx_rate * 100;
Percent_influx_totalT3_RB_blood_exit_rate = totalT3_RB_blood_exit_rate/RBT3_influx_rate  * 100;

% T4 tissue efflux rate (pmol/S)
LiverT4_efflux_rate = param.k30*T4LT*param.fuT4LT;
Percent_influx_LiverT4_efflux_rate = LiverT4_efflux_rate/LiverT4_influx_rate * 100;
RBT4_efflux_rate = param.k28*T4RBT*param.fuT4RBT;
Percent_influx_RBT4_efflux_rate = RBT4_efflux_rate/RBT4_influx_rate * 100;
BodyT4_efflux_rate = LiverT4_efflux_rate + RBT4_efflux_rate;

% T4 clearance rate (pmol/S)
LiverT4_clearance_rate = (param.k26+param.k34)*T4LT*param.fuT4LT*param.VLT;
Percent_influx_LiverT4_clearance_rate = LiverT4_clearance_rate/LiverT4_influx_rate * 100;
RBT4_clearance_rate = (param.k24+param.k32)*T4RBT*param.fuT4RBT*param.VRBT;
Percent_influx_RBT4_clearance_rate = RBT4_clearance_rate/RBT4_influx_rate * 100;
BodyT4_clearance_rate = LiverT4_clearance_rate + RBT4_clearance_rate;

Percent_LiverT4_clearance = LiverT4_clearance_rate / param.k20;
Percent_RBT4_clearance = RBT4_clearance_rate / param.k20;


% T3 tissue efflux rate (pmol/S)
LiverT3_efflux_rate = param.k31*T3LT*param.fuT3LT;
Percent_influx_LiverT3_efflux_rate = LiverT3_efflux_rate/LiverT3_influx_rate * 100;
RBT3_efflux_rate = param.k29*T3RBT*param.fuT3RBT;
Percent_influx_RBT3_efflux_rate = RBT3_efflux_rate/RBT3_influx_rate * 100;
BodyT3_efflux_rate = LiverT3_efflux_rate + RBT3_efflux_rate;

% T3 clearance rate
LiverT3_clearance_rate =param.k35*T3LT*param.fuT3LT*param.VLT;
Percent_influx_LiverT3_clearance_rate = LiverT3_clearance_rate/LiverT3_influx_rate * 100;
RBT3_clearance_rate = param.k33*T3RBT*param.fuT3RBT*param.VRBT;
Percent_influx_RBT3_clearance_rate = RBT3_clearance_rate/RBT3_influx_rate * 100;
BodyT3_clearance_rate = LiverT3_clearance_rate + RBT3_clearance_rate;

% T3 tissue production rate
LiverT3_production_rate = param.k26*T4LT*param.fuT4LT*param.VLT;
Percent_influx_LiverT3_production_rate = LiverT3_production_rate/LiverT3_influx_rate * 100;
RBT3_production_rate = param.k24*T4RBT*param.fuT4RBT*param.VRBT;
Percent_influx_RBT3_production_rate = RBT3_production_rate/RBT3_influx_rate * 100;
T3_production_rate = LiverT3_production_rate + RBT3_production_rate + param.k22;
Percent_LiverT3_production_rate = LiverT3_production_rate/T3_production_rate;
Percent_RBT3_production_rate = RBT3_production_rate/T3_production_rate;
Percent_ThyroidT3_production_rate = param.k22/T3_production_rate;

% T3 percentage clearance rate
Percent_LiverT3_clearance = LiverT3_clearance_rate / T3_production_rate;
Percent_RBT3_clearance = RBT3_clearance_rate / T3_production_rate;

% T4+T3 clearance rate
BodyT4T3_clearance_rate = param.k32*T4RBT*param.fuT4RBT*param.VRBT + param.k34*T4LT*param.fuT4LT*param.VLT + param.k33*T3RBT*param.fuT3RBT*param.VRBT + param.k35*T3LT*param.fuT3LT*param.VLT;

% T3 production rate in tissues
TissueT3_production_rate = (param.k24*T4RBT*param.fuT4RBT*param.VRBT + param.k26*T4LT*param.fuT4LT*param.VLT);

% Blood transit time through tissues
Transit_time_liver = param.VLB/param.QL;
Transit_time_RB = param.VRBB/param.QRB;

% T4 and T3 differnetial between arterial and venous blood concentrations
CACVfT4T = freeT4T - freeT4B;
CACVfT4L = freeT4L - freeT4B;
CACVfT4RB = freeT4RB - freeT4B;

CACVfT4T_percent = (freeT4T - freeT4B)/freeT4B * 100;
CACVfT4L_percent = (freeT4L - freeT4B)/freeT4B * 100;
CACVfT4RB_percent = (freeT4RB - freeT4B)/freeT4B * 100;

CACVfT3T = freeT3T - freeT3B;
CACVfT3L = freeT3L - freeT3B;
CACVfT3RB = freeT3RB - freeT3B;

CACVfT3T_percent = (freeT3T - freeT3B)/freeT3B * 100;
CACVfT3L_percent = (freeT3L - freeT3B)/freeT3B * 100;
CACVfT3RB_percent = (freeT3RB - freeT3B)/freeT3B * 100;

CACVT4TBGT = T4TBGT - T4TBGB;
CACVT4TBGL = T4TBGL - T4TBGB;
CACVT4TBGRB = T4TBGRB - T4TBGB;

CACVT4TBGT_percent = (T4TBGT - T4TBGB)/T4TBGB * 100;
CACVT4TBGL_percent = (T4TBGL - T4TBGB)/T4TBGB * 100;
CACVT4TBGRB_percent = (T4TBGRB - T4TBGB)/T4TBGB * 100;

CACVT3TBGT = T3TBGT - T3TBGB;
CACVT3TBGL = T3TBGL - T3TBGB;
CACVT3TBGRB = T3TBGRB - T3TBGB;

CACVT3TBGT_percent = (T3TBGT - T3TBGB)/T3TBGB * 100;
CACVT3TBGL_percent = (T3TBGL - T3TBGB)/T3TBGB * 100;
CACVT3TBGRB_percent = (T3TBGRB - T3TBGB)/T3TBGB * 100;

CACVT4TTRT = T4TTRT - T4TTRB;
CACVT4TTRL = T4TTRL - T4TTRB;
CACVT4TTRRB = T4TTRRB - T4TTRB;

CACVT4TTRT_percent = (T4TTRT - T4TTRB)/T4TTRB * 100;
CACVT4TTRL_percent = (T4TTRL - T4TTRB)/T4TTRB * 100;
CACVT4TTRRB_percent = (T4TTRRB - T4TTRB)/T4TTRB * 100;

CACVT3TTRT = T3TTRT - T3TTRB;
CACVT3TTRL = T3TTRL - T3TTRB;
CACVT3TTRRB = T3TTRRB - T3TTRB;

CACVT3TTRT_percent = (T3TTRT - T3TTRB)/T3TTRB * 100;
CACVT3TTRL_percent = (T3TTRL - T3TTRB)/T3TTRB * 100;
CACVT3TTRRB_percent = (T3TTRRB - T3TTRB)/T3TTRB * 100;

CACVT4ALBT = T4ALBT - T4ALBB;
CACVT4ALBL = T4ALBL - T4ALBB;
CACVT4ALBRB = T4ALBRB - T4ALBB;

CACVT4ALBT_percent = (T4ALBT - T4ALBB)/T4ALBB * 100;
CACVT4ALBL_percent = (T4ALBL - T4ALBB)/T4ALBB * 100;
CACVT4ALBRB_percent = (T4ALBRB - T4ALBB)/T4ALBB * 100;

CACVT3ALBT = T3ALBT - T3ALBB;
CACVT3ALBL = T3ALBL - T3ALBB;
CACVT3ALBRB = T3ALBRB - T3ALBB;

CACVT3ALBT_percent = (T3ALBT - T3ALBB)/T3ALBB * 100;
CACVT3ALBL_percent = (T3ALBL - T3ALBB)/T3ALBB * 100;
CACVT3ALBRB_percent = (T3ALBRB - T3ALBB)/T3ALBB * 100;

% Binding Rates in Body Blood
body_blood_T4TBG_assoc_rate = param.k1*freeT4B*TBGB;
body_blood_T4TBG_disassoc_rate = param.k2*T4TBGB;
body_blood_difference_T4TBG_rate = body_blood_T4TBG_disassoc_rate - body_blood_T4TBG_assoc_rate; 
body_blood_percent_difference_T4TBG_rate = body_blood_difference_T4TBG_rate/body_blood_T4TBG_assoc_rate*100;

body_blood_T4TTR_assoc_rate = param.k3*freeT4B*TTRB;
body_blood_T4TTR_disassoc_rate = param.k4*T4TTRB;
body_blood_difference_T4TTR_rate = body_blood_T4TTR_disassoc_rate - body_blood_T4TTR_assoc_rate;
body_blood_percent_difference_T4TTR_rate = body_blood_difference_T4TTR_rate/body_blood_T4TTR_assoc_rate*100;

body_blood_T4ALB_assoc_rate = param.k5*freeT4B*ALBB;
body_blood_T4ALB_disassoc_rate = param.k6*T4ALBB;
body_blood_difference_T4ALB_rate = body_blood_T4ALB_disassoc_rate - body_blood_T4ALB_assoc_rate;
body_blood_percent_difference_T4ALB_rate = body_blood_difference_T4ALB_rate/body_blood_T4ALB_assoc_rate*100;

body_blood_T3TBG_assoc_rate = param.k7*freeT3B*TBGB;
body_blood_T3TBG_disassoc_rate = param.k8*T3TBGB;
body_blood_difference_T3TBG_rate = body_blood_T3TBG_disassoc_rate - body_blood_T3TBG_assoc_rate; 
body_blood_percent_difference_T3TBG_rate = body_blood_difference_T3TBG_rate/body_blood_T3TBG_assoc_rate*100;

body_blood_T3TTR_assoc_rate = param.k9*freeT3B*TTRB;
body_blood_T3TTR_disassoc_rate = param.k10*T3TTRB;
body_blood_difference_T3TTR_rate = body_blood_T3TTR_disassoc_rate - body_blood_T3TTR_assoc_rate; 
body_blood_percent_difference_T3TTR_rate = body_blood_difference_T3TTR_rate/body_blood_T3TTR_assoc_rate*100;

body_blood_T3ALB_assoc_rate = param.k11*freeT3B*ALBB;
body_blood_T3ALB_disassoc_rate = param.k12*T3ALBB;
body_blood_difference_T3ALB_rate = body_blood_T3ALB_disassoc_rate - body_blood_T3ALB_assoc_rate; 
body_blood_percent_difference_T3ALB_rate = body_blood_difference_T3ALB_rate/body_blood_T3ALB_assoc_rate*100;

% Binding Rates in Thyroid blood
thyroid_T4TBG_assoc_rate = param.k1*freeT4T*TBGT;
thyroid_T4TBG_disassoc_rate = param.k2*T4TBGT;
thyroid_difference_T4TBG_rate = thyroid_T4TBG_disassoc_rate - thyroid_T4TBG_assoc_rate; 
thyroid_percent_difference_T4TBG_rate = thyroid_difference_T4TBG_rate/thyroid_T4TBG_assoc_rate*100;

thyroid_T4TTR_assoc_rate = param.k3*freeT4T*TTRT;
thyroid_T4TTR_disassoc_rate = param.k4*T4TTRT;
thyroid_difference_T4TTR_rate = thyroid_T4TTR_disassoc_rate - thyroid_T4TTR_assoc_rate; 
thyroid_percent_difference_T4TTR_rate = thyroid_difference_T4TTR_rate/thyroid_T4TTR_assoc_rate*100;

thyroid_T4ALB_assoc_rate = param.k5*freeT4T*ALBT;
thyroid_T4ALB_disassoc_rate = param.k6*T4ALBT;
thyroid_difference_T4ALB_rate = thyroid_T4ALB_disassoc_rate - thyroid_T4ALB_assoc_rate; 
thyroid_percent_difference_T4ALB_rate = thyroid_difference_T4ALB_rate/thyroid_T4ALB_assoc_rate*100;

thyroid_T3TBG_assoc_rate = param.k7*freeT3T*TBGT;
thyroid_T3TBG_disassoc_rate = param.k8*T3TBGT;
thyroid_difference_T3TBG_rate = thyroid_T3TBG_disassoc_rate - thyroid_T3TBG_assoc_rate; 
thyroid_percent_difference_T3TBG_rate = thyroid_difference_T3TBG_rate/thyroid_T3TBG_assoc_rate*100;

thyroid_T3TTR_assoc_rate = param.k9*freeT3T*TTRT;
thyroid_T3TTR_disassoc_rate = param.k10*T3TTRT;
thyroid_difference_T3TTR_rate = thyroid_T3TTR_disassoc_rate - thyroid_T3TTR_assoc_rate; 
thyroid_percent_difference_T3TTR_rate = thyroid_difference_T3TTR_rate/thyroid_T3TTR_assoc_rate*100;

thyroid_T3ALB_assoc_rate = param.k11*freeT3T*ALBT;
thyroid_T3ALB_disassoc_rate = param.k12*T3ALBT;
thyroid_difference_T3ALB_rate = thyroid_T3ALB_disassoc_rate - thyroid_T3ALB_assoc_rate; 
thyroid_percent_difference_T3ALB_rate = thyroid_difference_T3ALB_rate/thyroid_T3ALB_assoc_rate*100;

% Binding Rates in Liver blood
liver_T4TBG_assoc_rate = param.k1*freeT4L*TBGL;
liver_T4TBG_disassoc_rate = param.k2*T4TBGL;
liver_difference_T4TBG_rate = liver_T4TBG_disassoc_rate - liver_T4TBG_assoc_rate; 
liver_percent_difference_T4TBG_rate = liver_difference_T4TBG_rate/liver_T4TBG_assoc_rate*100;

liver_T4TTR_assoc_rate = param.k3*freeT4L*TTRL;
liver_T4TTR_disassoc_rate = param.k4*T4TTRL;
liver_difference_T4TTR_rate = liver_T4TTR_disassoc_rate - liver_T4TTR_assoc_rate; 
liver_percent_difference_T4TTR_rate = liver_difference_T4TTR_rate/liver_T4TTR_assoc_rate*100;

liver_T4ALB_assoc_rate = param.k5*freeT4L*ALBL;
liver_T4ALB_disassoc_rate = param.k6*T4ALBL;
liver_difference_T4ALB_rate = liver_T4ALB_disassoc_rate - liver_T4ALB_assoc_rate; 
liver_percent_difference_T4ALB_rate = liver_difference_T4ALB_rate/liver_T4ALB_assoc_rate*100

liver_T3TBG_assoc_rate = param.k7*freeT3L*TBGL;
liver_T3TBG_disassoc_rate = param.k8*T3TBGL;
liver_difference_T3TBG_rate = liver_T3TBG_disassoc_rate - liver_T3TBG_assoc_rate; 
liver_percent_difference_T3TBG_rate = liver_difference_T3TBG_rate/liver_T3TBG_assoc_rate*100;

liver_T3TTR_assoc_rate = param.k9*freeT3L*TTRL;
liver_T3TTR_disassoc_rate = param.k10*T3TTRL;
liver_difference_T3TTR_rate = liver_T3TTR_disassoc_rate - liver_T3TTR_assoc_rate; 
liver_percent_difference_T3TTR_rate = liver_difference_T3TTR_rate/liver_T3TTR_assoc_rate*100;

liver_T3ALB_assoc_rate = param.k11*freeT3L*ALBL;
liver_T3ALB_disassoc_rate = param.k12*T3ALBL;
liver_difference_T3ALB_rate = liver_T3ALB_disassoc_rate - liver_T3ALB_assoc_rate;
liver_percent_difference_T3ALB_rate = liver_difference_T3ALB_rate/liver_T3ALB_assoc_rate*100;

% Binding Rates in RB blood
RB_T4TBG_assoc_rate = param.k1*freeT4RB*TBGRB;
RB_T4TBG_disassoc_rate = param.k2*T4TBGRB;
RB_difference_T4TBG_rate = RB_T4TBG_disassoc_rate - RB_T4TBG_assoc_rate; 
RB_percent_difference_T4TBG_rate = RB_difference_T4TBG_rate/RB_T4TBG_assoc_rate*100;

RB_T4TTR_assoc_rate = param.k3*freeT4RB*TTRRB;
RB_T4TTR_disassoc_rate = param.k4*T4TTRRB;
RB_difference_T4TTR_rate = RB_T4TTR_disassoc_rate - RB_T4TTR_assoc_rate; 
RB_percent_difference_T4TTR_rate = RB_difference_T4TTR_rate/RB_T4TTR_assoc_rate*100;

RB_T4ALB_assoc_rate = param.k5*freeT4RB*ALBRB;
RB_T4ALB_disassoc_rate = param.k6*T4ALBRB;
RB_difference_T4ALB_rate = RB_T4ALB_disassoc_rate - RB_T4ALB_assoc_rate; 
RB_percent_difference_T4ALB_rate = RB_difference_T4ALB_rate/RB_T4ALB_assoc_rate*100;

RB_T3TBG_assoc_rate = param.k7*freeT3RB*TBGRB;
RB_T3TBG_disassoc_rate = param.k8*T3TBGRB;
RB_difference_T3TBG_rate = RB_T3TBG_disassoc_rate - RB_T3TBG_assoc_rate; 
RB_percent_difference_T3TBG_rate = RB_difference_T3TBG_rate/RB_T3TBG_assoc_rate*100;

RB_T3TTR_assoc_rate = param.k9*freeT3RB*TTRRB;
RB_T3TTR_disassoc_rate = param.k10*T3TTRRB;
RB_difference_T3TTR_rate = RB_T3TTR_disassoc_rate - RB_T3TTR_assoc_rate; 
RB_percent_difference_T3TTR_rate = RB_difference_T3TTR_rate/RB_T3TTR_assoc_rate*100;

RB_T3ALB_assoc_rate = param.k11*freeT3RB*ALBRB;
RB_T3ALB_disassoc_rate = param.k12*T3ALBRB;
RB_difference_T3ALB_rate = RB_T3ALB_disassoc_rate - RB_T3ALB_assoc_rate; 
RB_percent_difference_T3ALB_rate = RB_difference_T3ALB_rate/RB_T3ALB_assoc_rate*100;

% 


%% -------------------------- Mass-balance check ---------------------------------%%

% T4 and T3 amounts in Entire Body at steady state
T4_in_entire_body_at_steady_state = totalT4_in_entire_body;
T3_in_entire_body_at_steady_state = totalT3_in_entire_body;

% T4 mass-balance check: by setting k20=0
y0 = y(end, :); %Using the steady-state values from the "Run to Steady State" section above as the initial values
param = default_param;
param.k20 = 0;
sampling_interval = 10;
tspan0 = [0:sampling_interval:60*24*3600]; %Running for 60 days
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t1,y1] = ode15s('TH_PBK_Nonspatial_ODE',tspan0, y0, options, param);

day = t1/3600/24;

%Simulation results
model.fT4B    = y1(:,1);
model.T4TBGB  = y1(:,2);
model.T4TTRB  = y1(:,3);
model.T4ALBB  = y1(:,4);

model.fT4T    = y1(:,12);
model.T4TBGT  = y1(:,13);
model.T4TTRT  = y1(:,14);
model.T4ALBT  = y1(:,15);

model.fT4RB    = y1(:,23);  
model.T4TBGRB  = y1(:,24);
model.T4TTRRB  = y1(:,25);
model.T4ALBRB  = y1(:,26);
model.T4RBT    = y1(:,34);

model.fT4L    = y1(:,36);
model.T4TBGL  = y1(:,37);
model.T4TTRL  = y1(:,38);
model.T4ALBL  = y1(:,39);
model.T4LT    = y1(:,44);

% Total T4 amounts in blood compartments
totalT4_in_Body_Blood = (model.fT4B + model.T4TBGB + model.T4TTRB + model.T4ALBB) * param.VB;
totalT4_in_Thyroid_blood = (model.fT4T + model.T4TBGT + model.T4TTRT + model.T4ALBT) * param.VTB;
totalT4_in_RB_blood = (model.fT4RB + model.T4TBGRB + model.T4TTRRB + model.T4ALBRB) * param.VRBB;
totalT4_in_Liver_blood = (model.fT4L + model.T4TBGL + model.T4TTRL + model.T4ALBL) * param.VLB;

% T4 amounts in tissue compartments
T4_in_RB_tissue = model.T4RBT * param.VRBT;
T4_in_Liver_tissue = model.T4LT * param.VLT;

% Cumulative T4 amounts metabolized in tissue compartments
T4_metabolized_in_RB_tissue = [];
T4_metabolized_in_Liver_tissue = [];
T4_metabolized_in_RB_tissue(1) = 0;
T4_metabolized_in_Liver_tissue(1) = 0;
for i=2:length(t1)
    T4_metabolized_in_RB_tissue(i) = T4_metabolized_in_RB_tissue(i-1) + model.T4RBT(i-1) * param.fuT4RBT * (param.k24+param.k32) * param.VRBT * sampling_interval;
    T4_metabolized_in_Liver_tissue(i) = T4_metabolized_in_Liver_tissue(i-1) + model.T4LT(i-1) * param.fuT4LT * (param.k26+param.k34) * param.VLT * sampling_interval;
end

% Sum of T4 amounts in all compartments and metabolized
T4_summed = totalT4_in_Body_Blood + totalT4_in_Thyroid_blood + totalT4_in_RB_blood + totalT4_in_Liver_blood + T4_in_RB_tissue + T4_in_Liver_tissue + T4_metabolized_in_RB_tissue' + T4_metabolized_in_Liver_tissue';

%Plotting 
figure(100)
hold on
plot([0,day(end)],[T4_in_entire_body_at_steady_state, T4_in_entire_body_at_steady_state])
plot(day(1:1000:end),totalT4_in_Body_Blood(1:1000:end))
plot(day(1:1000:end),totalT4_in_Thyroid_blood(1:1000:end))
plot(day(1:1000:end),totalT4_in_RB_blood(1:1000:end))
plot(day(1:1000:end),totalT4_in_Liver_blood(1:1000:end))
plot(day(1:1000:end),T4_in_RB_tissue(1:1000:end))
plot(day(1:1000:end),T4_in_Liver_tissue(1:1000:end))
plot(day(1:1000:end),T4_metabolized_in_RB_tissue(1:1000:end))
plot(day(1:1000:end),T4_metabolized_in_Liver_tissue(1:1000:end))
plot(day(1:1000:end),T4_summed(1:1000:end))
xlabel('Day')
ylabel('T4 amount (pmol)')
box on
pbaspect([2,1,1])


% T3 mass-balance check: by setting k22=0, k24=0, and k26=0
y0 = y(end, :); %Using the steady-state values from the "Run to Steady State" section above as the initial values
param = default_param;
param.k22 = 0;
param.k32 = param.k32 + param.k24; %To keep RB T4 metabolism the same
param.k24 = 0;
param.k34 = param.k34 + param.k26; %To keep Liver T4 metabolism the same
param.k26 = 0;
sampling_interval = 10;
tspan0 = [0:sampling_interval:10*24*3600]; %Running for 10 days
options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t2,y2] = ode15s('TH_PBK_Nonspatial_ODE',tspan0, y0, options, param);

day = t2/3600/24;

%Simulation results
model.fT3B    = y2(:,5);
model.T3TBGB  = y2(:,6);
model.T3TTRB  = y2(:,7);
model.T3ALBB  = y2(:,8);

model.fT3T    = y2(:,16);
model.T3TBGT  = y2(:,17);
model.T3TTRT  = y2(:,18);
model.T3ALBT  = y2(:,19);

model.fT3RB    = y2(:,27);  
model.T3TBGRB  = y2(:,28);
model.T3TTRRB  = y2(:,29);
model.T3ALBRB  = y2(:,30);
model.T3RBT    = y2(:,35);

model.fT3L    = y2(:,40);
model.T3TBGL  = y2(:,41);
model.T3TTRL  = y2(:,42);
model.T3ALBL  = y2(:,43);
model.T3LT    = y2(:,45);

% Total T3 amounts in blood compartments
totalT3_in_Body_Blood = (model.fT3B + model.T3TBGB + model.T3TTRB + model.T3ALBB) * param.VB;
totalT3_in_Thyroid_blood = (model.fT3T + model.T3TBGT + model.T3TTRT + model.T3ALBT) * param.VTB;
totalT3_in_RB_blood = (model.fT3RB + model.T3TBGRB + model.T3TTRRB + model.T3ALBRB) * param.VRBB;
totalT3_in_Liver_blood = (model.fT3L + model.T3TBGL + model.T3TTRL + model.T3ALBL) * param.VLB;

% T3 amounts in tissue compartments
T3_in_RB_tissue = model.T3RBT * param.VRBT;
T3_in_Liver_tissue = model.T3LT * param.VLT;

% Cumulative T3 amounts metabolized in tissue compartments
T3_metabolized_in_RB_tissue = [];
T3_metabolized_in_Liver_tissue = [];
T3_metabolized_in_RB_tissue(1) = 0;
T3_metabolized_in_Liver_tissue(1) = 0;
for i=2:length(t2)
    T3_metabolized_in_RB_tissue(i) = T3_metabolized_in_RB_tissue(i-1) + model.T3RBT(i-1) * param.fuT3RBT * param.k33 * param.VRBT * sampling_interval;
    T3_metabolized_in_Liver_tissue(i) = T3_metabolized_in_Liver_tissue(i-1) + model.T3LT(i-1) * param.fuT3LT * param.k35 * param.VLT * sampling_interval;
end

% Sum of T3 amounts in all compartments and metabolized
T3_summed = totalT3_in_Body_Blood + totalT3_in_Thyroid_blood + totalT3_in_RB_blood + totalT3_in_Liver_blood + T3_in_RB_tissue + T3_in_Liver_tissue + T3_metabolized_in_RB_tissue' + T3_metabolized_in_Liver_tissue';

%Plotting 
figure(200)
hold on
plot([0,day(end)],[T3_in_entire_body_at_steady_state, T3_in_entire_body_at_steady_state])
plot(day(1:1000:end),totalT3_in_Body_Blood(1:1000:end))
plot(day(1:1000:end),totalT3_in_Thyroid_blood(1:1000:end))
plot(day(1:1000:end),totalT3_in_RB_blood(1:1000:end))
plot(day(1:1000:end),totalT3_in_Liver_blood(1:1000:end))
plot(day(1:1000:end),T3_in_RB_tissue(1:1000:end))
plot(day(1:1000:end),T3_in_Liver_tissue(1:1000:end))
plot(day(1:1000:end),T3_metabolized_in_RB_tissue(1:1000:end))
plot(day(1:1000:end),T3_metabolized_in_Liver_tissue(1:1000:end))
plot(day(1:1000:end),T3_summed(1:1000:end))
xlabel('Day')
ylabel('T3 amount (pmol)')
box on
pbaspect([2,1,1])


%% -------------------------- Sensitivity analysis (Table 6) ---------------------------------%%
y0 = y(end, :); %Using the steady-state values from the "Run to Steady State" section above as the initial values
param = default_param;

%Obtain basal steady-state values
ss_freeT4 = y0(1);
ss_freeT3 = y0(5);
ss_T4RBT = y0(34);
ss_T3RBT = y0(35);
ss_T4LT = y0(44);
ss_T3LT = y0(45);

tspan1 = [0:1000:1000*24*3600]; %Running for 1000 days

percent_change = 0.01;

%Increase the parameter value
    param.k24 = default_param.k24 * (1+percent_change); %replace k24 here for other parameters: k26,k32,k33,k34,k35
    options = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [t1,y1] = ode15s('TH_PBK_Nonspatial_ODE_clamped',tspan1, y0, options, param);

    %Obtain steady-state values after parameter change
    positive_ss_freeT4  = y1(end,1);
    positive_ss_freeT3  = y1(end,5);
    positive_ss_T4RBT   = y1(end,34);
    positive_ss_T3RBT   = y1(end,35);
    positive_ss_T4LT    = y1(end,44);
    positive_ss_T3LT    = y1(end,45);

%Decrease the parameter value
    param.k24 = default_param.k24 * (1-percent_change); %replace k24 here for other parameters: k26,k32,k33,k34,k35
    options = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [t1,y1] = ode15s('TH_PBK_Nonspatial_ODE_clamped',tspan1, y0, options, param);

    %Obtain steady-state values after parameter change
    negative_ss_freeT4  = y1(end,1);
    negative_ss_freeT3  = y1(end,5);
    negative_ss_T4RBT   = y1(end,34);
    negative_ss_T3RBT   = y1(end,35);
    negative_ss_T4LT    = y1(end,44);
    negative_ss_T3LT    = y1(end,45);

%Calculate sensitivity coefficient
    sensitivity_coeff_freeT4 = mean([(positive_ss_freeT4 - ss_freeT4)/ss_freeT4/percent_change, (ss_freeT4 - negative_ss_freeT4)/ss_freeT4/percent_change])
    sensitivity_coeff_freeT3 = mean([(positive_ss_freeT3 - ss_freeT3)/ss_freeT3/percent_change, (ss_freeT3 - negative_ss_freeT3)/ss_freeT3/percent_change])
    sensitivity_coeff_T4RBT = mean([(positive_ss_T4RBT - ss_T4RBT)/ss_T4RBT/percent_change, (ss_T4RBT - negative_ss_T4RBT)/ss_T4RBT/percent_change])
    sensitivity_coeff_T3RBT = mean([(positive_ss_T3RBT - ss_T3RBT)/ss_T3RBT/percent_change, (ss_T3RBT - negative_ss_T3RBT)/ss_T3RBT/percent_change])
    sensitivity_coeff_T4LT = mean([(positive_ss_T4LT - ss_T4LT)/ss_T4LT/percent_change, (ss_T4LT - negative_ss_T4LT)/ss_T4LT/percent_change])
    sensitivity_coeff_T3LT = mean([(positive_ss_T3LT - ss_T3LT)/ss_T3LT/percent_change, (ss_T3LT - negative_ss_T3LT)/ss_T3LT/percent_change])


%% -------------- Steady-state dose-response to continuouse X exposure (Fig. 10)-------------%%
y0 = y(end, :); %Using the steady-state values from the section above as the initial values
param = default_param;

X_total_halflife = 1; %half-life of total X in hours

% To turn on X binding to TBG (Fig. 10A-10C)
param.k43 = 0.018;
param.k42 = param.k43/param.kdXTBG;
X_halflife = X_total_halflife/3894*param.VB/VPtot; %X_halflife refers to half-life of free X in Body Blood (VB) in hours assuming X is only produced and metabolized in Body Blood; 3894 is the ratio of XTBG to free X when free X = 15 pM. 

% To turn on X binding to TTR (Fig. 10D-10F)
% param.k39 = 0.0832;
% param.k38 = param.k39/param.kdXTTR;
% X_halflife = X_total_halflife/1064*param.VB/VPtot; %X_halflife refers to half-life of free X in Body Blood (VB) in hours assuming X is only produced and metabolized in Body Blood; 1064 is the ratio of XTTR to free X when free X = 15 pM. 

XBs = [];
fT4s = [];
fT3s = [];
T4TBGBs = [];
T4TTRBs = [];
T4ALBBs = [];
T3TBGBs = [];
T3TTRBs = [];
T3ALBBs = [];
XTBGBs = [];
XTTRBs = [];
totalT4Bs = [];
totalT3Bs = [];

percent_increase = 1.5;
i=0.01;

while i<=1.5E6
    i
    %Parameters for X production and clearance in the Body Blood compartment(VB)
    param.k37 = log(2)/(X_halflife*3600);
    param.k36 = i * param.VB * param.k37;
    
    tspan1 = [0:1000:100*24*3600]; %Running for 100 days

    options = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [t1,y1] = ode15s('TH_PBK_Nonspatial_ODE',tspan1, y0, options, param);

%Simulation result
    model.fT4B    = y1(:,1);
    model.T4TBGB  = y1(:,2);
    model.T4TTRB  = y1(:,3);
    model.T4ALBB  = y1(:,4);

    model.fT3B    = y1(:,5);
    model.T3TBGB  = y1(:,6);
    model.T3TTRB  = y1(:,7);
    model.T3ALBB  = y1(:,8);

    model.TBGB    = y1(:,9);
    model.TTRB    = y1(:,10);
    model.ALBB    = y1(:,11);

    model.fT4T    = y1(:,12);
    model.T4TBGT  = y1(:,13);
    model.T4TTRT  = y1(:,14);
    model.T4ALBT  = y1(:,15);
    model.fT3T    = y1(:,16);
    model.T3TBGT  = y1(:,17);
    model.T3TTRT  = y1(:,18);
    model.T3ALBT  = y1(:,19);
    model.TBGT    = y1(:,20);
    model.TTRT    = y1(:,21);
    model.ALBT    = y1(:,22);

    model.fT4RB    = y1(:,23);  
    model.T4TBGRB  = y1(:,24);
    model.T4TTRRB  = y1(:,25);
    model.T4ALBRB  = y1(:,26);
    model.fT3RB    = y1(:,27);
    model.T3TBGRB  = y1(:,28);
    model.T3TTRRB  = y1(:,29);
    model.T3ALBRB  = y1(:,30);
    model.TBGRB    = y1(:,31);
    model.TTRRB    = y1(:,32);
    model.ALBRB    = y1(:,33);
    model.T4RBT    = y1(:,34);
    model.T3RBT    = y1(:,35);

    model.fT4L    = y1(:,36);
    model.T4TBGL  = y1(:,37);
    model.T4TTRL  = y1(:,38);
    model.T4ALBL  = y1(:,39);
    model.fT3L    = y1(:,40);
    model.T3TBGL  = y1(:,41);
    model.T3TTRL  = y1(:,42);
    model.T3ALBL  = y1(:,43);
    model.T4LT    = y1(:,44);
    model.T3LT    = y1(:,45);

    model.XB      = y1(:,46);
    model.XTTRB   = y1(:,47);
    model.XT      = y1(:,48);
    model.XTTRT   = y1(:,49);
    model.XRB     = y1(:,50);
    model.XTTRRB  = y1(:,51);
    model.XL      = y1(:,52);
    model.XTTRL   = y1(:,53);
    model.XTBGB   = y1(:,54);
    model.XTBGT   = y1(:,55);
    model.XTBGRB  = y1(:,56);
    model.XTBGL   = y1(:,57);

    model.TBGL = (param.TBGtot - (model.TBGB + model.T4TBGB + model.T3TBGB + model.XTBGB)*param.VB - (model.TBGT + model.T4TBGT + model.T3TBGT + model.XTBGT)*param.VTB - (model.TBGRB + model.T4TBGRB + model.T3TBGRB + model.XTBGRB)*param.VRBB - (model.T4TBGL + model.T3TBGL + model.XTBGL)*param.VLB)/param.VLB;
    model.TTRL = (param.TTRtot - (model.TTRB + model.T4TTRB + model.T3TTRB + model.XTTRB)*param.VB - (model.TTRT + model.T4TTRT + model.T3TTRT + model.XTTRT)*param.VTB - (model.TTRRB + model.T4TTRRB + model.T3TTRRB + model.XTTRRB)*param.VRBB - (model.T4TTRL + model.T3TTRL + model.XTTRL)*param.VLB)/param.VLB;
    model.ALBL = (param.ALBtot - (model.ALBB + model.T4ALBB + model.T3ALBB)*param.VB - (model.ALBT + model.T4ALBT + model.T3ALBT)*param.VTB - (model.ALBRB + model.T4ALBRB + model.T3ALBRB)*param.VRBB - (model.T4ALBL + model.T3ALBL)*param.VLB)/param.VLB;
    
    freeT4B = model.fT4B(end);
    T4TBGB = model.T4TBGB(end);
    T4TTRB = model.T4TTRB(end);
    T4ALBB = model.T4ALBB(end);
    freeT3B = model.fT3B(end);
    T3TBGB = model.T3TBGB(end);
    T3TTRB = model.T3TTRB(end);
    T3ALBB = model.T3ALBB(end);
    XTBGB  = model.XTBGB(end);
    XTTRB  = model.XTTRB(end);
    
    boundT4B = T4TBGB + T4TTRB + T4ALBB;
    totalT4B = boundT4B + freeT4B;
    
    boundT3B = T3TBGB + T3TTRB + T3ALBB;
    totalT3B = boundT3B + freeT3B;
    
    XBs = [XBs, model.XB(end)];
    fT4s = [fT4s, freeT4B];
    fT3s = [fT3s, freeT3B];
    T4TBGBs = [T4TBGBs, T4TBGB];
    T4TTRBs = [T4TTRBs, T4TTRB];
    T4ALBBs = [T4ALBBs, T4ALBB];
    T3TBGBs = [T3TBGBs, T3TBGB];
    T3TTRBs = [T3TTRBs, T3TTRB];
    T3ALBBs = [T3ALBBs, T3ALBB];
    XTBGBs  = [XTBGBs, XTBGB];
    XTTRBs  = [XTTRBs, XTTRB];
    totalT4Bs = [totalT4Bs, totalT4B];
    totalT3Bs = [totalT3Bs, totalT3B];
    
    
    i = i*percent_increase;
end


% Fig. 10A: XTBG vs. free X curve
figure(101)
semilogx(XBs, XTBGBs);
hold on
xlabel('Free X')
ylabel('XTBG')

% Fig. 10B: T4TBG and Total T4 vs. free X curve
figure(102)
semilogx(XBs, T4TBGBs);
hold on
semilogx(XBs, totalT4Bs);
xlabel('Free X')
ylabel('T4TBG or Total T4')

% Fig. 10C: T3TBG and Total T3 vs. free X curve
figure(103)
semilogx(XBs, T3TBGBs);
hold on
semilogx(XBs, totalT3Bs);
xlabel('Free X')
ylabel('T3TBG or Total T3')

% Fig. 10D: XTTR vs. free X curve
figure(104)
semilogx(XBs, XTTRBs);
hold on
xlabel('Free X')
ylabel('XTTR')

% Fig. 10E: T4TTR and Total T4 vs. free X curve
figure(105)
semilogx(XBs, T4TTRBs);
hold on
semilogx(XBs, totalT4Bs);
xlabel('Free X')
ylabel('T4TTR or Total T4')

% Fig. 10F: T3TTR and Total T3 vs. free X curve
figure(106)
semilogx(XBs, T3TTRBs);
hold on
semilogx(XBs, totalT3Bs);
xlabel('Free X')
ylabel('T3TTR or Total T3 ')


%% -------------- TH response to intermittent daily X exposure (Figs. 11, S16-S20)--------------%%
param = default_param;

exposure_time = 8; %Daily exposure time in hours. 

X_total_halflife = 1; % Change it to 1 for Figs. 11 & S18, to 10 for Figs. S16 & S19, to 100 for Figs. S17 & S20


% To turn on X binding to TBG (Figs. 11, S16-S17)
param.k43 = 0.018;
param.k42 = param.k43/param.kdXTBG;
X_halflife = X_total_halflife/3894*param.VB/VPtot; %X_halflife refers to half-life of free X in Body Blood (VB) in hours assuming X is only produced and metabolized in Body Blood; 3894 is the ratio of XTBG to free X when free X = 15 pM. 

% To turn on X binding to TTR (Figs. S18-S20)
% param.k39 = 0.0832;
% param.k38 = param.k39/param.kdXTTR;
% X_halflife = X_total_halflife/1064*param.VB/VPtot; %X_halflife refers to half-life of free X in Body Blood (VB) in hours assuming X is only produced and metabolized in Body Blood; 1064 is the ratio of XTTR to free X when free X = 15 pM. 

param.k37 = log(2)/(X_halflife*3600);

y_output = [];
range = [1:1:60]; % Number of days of exposure to X
ts = [];
for x = range
    %Turning on exposure to X
    if x == 1
        y0 = y(end, :); %On day 1, use the steady-state values from the "Run to Steady State" section
    else
        y0 = y_output(end,:); %On remaining days, use the last time points on the previous day
    end
    tspan1 = [(x-1)*24*3600 : 100 : ((x-1)*24 + exposure_time)*3600];
    param.k36 = 15*param.VB*param.k37; %Turning on exposure to X
    options = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [t1,y1] = ode15s('TH_PBK_Nonspatial_ODE',tspan1, y0, options, param);
    
    if x == 1
        ts = [ts; t1];
        y_output = [y_output; y1];
    else
        ts = [ts; t1(2:end)];
        y_output = [y_output; y1(2:end,:)];
    end

    %Turning off exposure to X
    y0 = y1(end, :);
    tspan2 = [((x-1)*24 + exposure_time)*3600 : 100 : x*24*3600];
    options = odeset('RelTol',1e-8,'AbsTol',1e-10);
    param.k36 = 0; %Turning off exposure to X
    [t2,y2] = ode15s('TH_PBK_Nonspatial_ODE',tspan2, y0, options, param);

    ts = [ts; t2(2:end)];
    y_output = [y_output; y2(2:end,:)];
 
end

day = ts/3600/24;
hour = ts/3600;

model.fT4B    = y_output(:,1);
model.T4TBGB  = y_output(:,2);
model.T4TTRB  = y_output(:,3);
model.T4ALBB  = y_output(:,4);

model.fT3B    = y_output(:,5);
model.T3TBGB  = y_output(:,6);
model.T3TTRB  = y_output(:,7);
model.T3ALBB  = y_output(:,8);

model.TBGB    = y_output(:,9);
model.TTRB    = y_output(:,10);
model.ALBB    = y_output(:,11);

model.fT4T    = y_output(:,12);
model.T4TBGT  = y_output(:,13);
model.T4TTRT  = y_output(:,14);
model.T4ALBT  = y_output(:,15);
model.fT3T    = y_output(:,16);
model.T3TBGT  = y_output(:,17);
model.T3TTRT  = y_output(:,18);
model.T3ALBT  = y_output(:,19);
model.TBGT    = y_output(:,20);
model.TTRT    = y_output(:,21);
model.ALBT    = y_output(:,22);

model.fT4RB    = y_output(:,23);  
model.T4TBGRB  = y_output(:,24);
model.T4TTRRB  = y_output(:,25);
model.T4ALBRB  = y_output(:,26);
model.fT3RB    = y_output(:,27);
model.T3TBGRB  = y_output(:,28);
model.T3TTRRB  = y_output(:,29);
model.T3ALBRB  = y_output(:,30);
model.TBGRB    = y_output(:,31);
model.TTRRB    = y_output(:,32);
model.ALBRB    = y_output(:,33);
model.T4RBT    = y_output(:,34);
model.T3RBT    = y_output(:,35);

model.fT4L    = y_output(:,36);
model.T4TBGL  = y_output(:,37);
model.T4TTRL  = y_output(:,38);
model.T4ALBL  = y_output(:,39);
model.fT3L    = y_output(:,40);
model.T3TBGL  = y_output(:,41);
model.T3TTRL  = y_output(:,42);
model.T3ALBL  = y_output(:,43);
model.T4LT    = y_output(:,44);
model.T3LT    = y_output(:,45);

model.XB      = y_output(:,46);
model.XTTRB   = y_output(:,47);
model.XT      = y_output(:,48);
model.XTTRT   = y_output(:,49);
model.XRB     = y_output(:,50);
model.XTTRRB  = y_output(:,51);
model.XL      = y_output(:,52);
model.XTTRL   = y_output(:,53);
model.XTBGB   = y_output(:,54);
model.XTBGT   = y_output(:,55);
model.XTBGRB  = y_output(:,56);
model.XTBGL   = y_output(:,57);

model.TBGL = (param.TBGtot - (model.TBGB + model.T4TBGB + model.T3TBGB + model.XTBGB)*param.VB - (model.TBGT + model.T4TBGT + model.T3TBGT + model.XTBGT)*param.VTB - (model.TBGRB + model.T4TBGRB + model.T3TBGRB + model.XTBGRB)*param.VRBB - (model.T4TBGL + model.T3TBGL + model.XTBGL)*param.VLB)/param.VLB;
model.TTRL = (param.TTRtot - (model.TTRB + model.T4TTRB + model.T3TTRB + model.XTTRB)*param.VB - (model.TTRT + model.T4TTRT + model.T3TTRT + model.XTTRT)*param.VTB - (model.TTRRB + model.T4TTRRB + model.T3TTRRB + model.XTTRRB)*param.VRBB - (model.T4TTRL + model.T3TTRL + model.XTTRL)*param.VLB)/param.VLB;
model.ALBL = (param.ALBtot - (model.ALBB + model.T4ALBB + model.T3ALBB)*param.VB - (model.ALBT + model.T4ALBT + model.T3ALBT)*param.VTB - (model.ALBRB + model.T4ALBRB + model.T3ALBRB)*param.VRBB - (model.T4ALBL + model.T3ALBL)*param.VLB)/param.VLB;


xlim_lower = 0;
xlim_upper = 40;
x_ratio = 5;
y_ratio = 1;
z_ratio = 1;
figure(110)
subplot_total = 8;
subplot(subplot_total,2,1);
plot(day, model.XB)

%For T4 related variables
xlim([xlim_lower xlim_upper])
pbaspect([x_ratio y_ratio z_ratio])
legend("X")
subplot(subplot_total,2,2);
plot(day, model.XTBGB)
xlim([xlim_lower xlim_upper])
pbaspect([x_ratio y_ratio z_ratio])
legend("XTBG")
subplot(subplot_total,2,3);
plot(day, model.fT4B)
xlim([xlim_lower xlim_upper])
pbaspect([x_ratio y_ratio z_ratio])
legend("fT4")
subplot(subplot_total,2,4);
plot(day, model.T4TBGB)
xlim([xlim_lower xlim_upper])
pbaspect([x_ratio y_ratio z_ratio])
legend("T4TBG")
subplot(subplot_total,2,5);
plot(day, model.T4TTRB);
xlim([xlim_lower xlim_upper])
pbaspect([x_ratio y_ratio z_ratio])
legend("T4TTR")
subplot(subplot_total,2,6);
plot(day, model.T4ALBB);
xlim([xlim_lower xlim_upper])
pbaspect([x_ratio y_ratio z_ratio])
legend("T4ALB")
subplot(subplot_total,2,7)
plot(day, model.fT4B+model.T4TBGB+model.T4TTRB+model.T4ALBB);
xlim([xlim_lower xlim_upper])
pbaspect([x_ratio y_ratio z_ratio])
legend("Total T4")
subplot(subplot_total,2,8);
plot(day, model.T4LT);
xlim([xlim_lower xlim_upper])
pbaspect([x_ratio y_ratio z_ratio])
legend("Liver T4")
subplot(subplot_total,2,9);
plot(day, model.T4RBT);
xlim([xlim_lower xlim_upper])
pbaspect([x_ratio y_ratio z_ratio])
legend("RB T4")

%For T3 related variables
subplot(subplot_total,2,10);
plot(day, model.fT3B)
xlim([xlim_lower xlim_upper])
pbaspect([x_ratio y_ratio z_ratio])
legend("fT3")
subplot(subplot_total,2,11);
plot(day, model.T3TBGB)
xlim([xlim_lower xlim_upper])
pbaspect([x_ratio y_ratio z_ratio])
legend("T3TBG")
subplot(subplot_total,2,12);
plot(day, model.T3TTRB);
xlim([xlim_lower xlim_upper])
pbaspect([x_ratio y_ratio z_ratio])
legend("T3TTR")
subplot(subplot_total,2,13);
plot(day, model.T3ALBB);
xlim([xlim_lower xlim_upper])
pbaspect([x_ratio y_ratio z_ratio])
legend("T3ALB")
subplot(subplot_total,2,14)
plot(day, model.fT3B+model.T3TBGB+model.T3TTRB+model.T3ALBB);
xlim([xlim_lower xlim_upper])
pbaspect([x_ratio y_ratio z_ratio])
legend("Total T3")
subplot(subplot_total,2,15);
plot(day, model.T3LT);
xlim([xlim_lower xlim_upper])
pbaspect([x_ratio y_ratio z_ratio])
legend("Liver T3")
xlabel("Day")
subplot(subplot_total,2,16);
plot(day, model.T3RBT);
xlim([xlim_lower xlim_upper])
pbaspect([x_ratio y_ratio z_ratio])
legend("RB T3")
xlabel("Day")


%% -------------- TH response to intermittent daily X exposure (Fig. 12)--------------%%
param = default_param;

X_total_halflife = 1;

% To turn on X binding to TBG 
param.k43 = 0.018;
param.k42 = param.k43/param.kdXTBG;
X_halflife = X_total_halflife/3894*param.VB/VPtot; %X_halflife refers to half-life of free X in Body Blood (VB) in hours assuming X is only produced and metabolized in Body Blood; 3894 is the ratio of XTBG to free X when free X = 15 pM. 

param.k37 = log(2)/(X_halflife*3600);

for i = [24,23,20,16,8,2]
    exposure_time = i %Daily exposure time in hours. 
    y_output = [];
    range = [1:1:52]; % Number of days of exposure to X
    ts = [];
    for x = range
        %Turning on exposure to X
        if x == 1
            y0 = y(end, :); %On day 1, use the steady-state values from the "Run to Steady State" section
        else
            y0 = y_output(end,:); %On remaining days, use the last time points on the previous day
        end
        tspan1 = [(x-1)*24*3600 : 100 : ((x-1)*24 + exposure_time)*3600];
        param.k36 = 15*param.VB*param.k37; %Turning on exposure to X
        options = odeset('RelTol',1e-8,'AbsTol',1e-10);
        [t1,y1] = ode15s('TH_PBK_Nonspatial_ODE',tspan1, y0, options, param);

        if x == 1
            ts = [ts; t1];
            y_output = [y_output; y1];
        else
            ts = [ts; t1(2:end)];
            y_output = [y_output; y1(2:end,:)];
        end

        %Turning off exposure to X
        if exposure_time ~= 24 %Avoid running simulation when the exposure duration is 24 h, which will produce an error
            y0 = y1(end, :);
            tspan2 = [((x-1)*24 + exposure_time)*3600 : 100 : x*24*3600];
            options = odeset('RelTol',1e-8,'AbsTol',1e-10);
            param.k36 = 0; %Turning off exposure to X
            [t2,y2] = ode15s('TH_PBK_Nonspatial_ODE',tspan2, y0, options, param);

            ts = [ts; t2(2:end)];
            y_output = [y_output; y2(2:end,:)];
        end

    end

    day = ts/3600/24;
    hour = ts/3600;

    model.fT4B    = y_output(:,1);
    model.T4TBGB  = y_output(:,2);
    model.T4TTRB  = y_output(:,3);
    model.T4ALBB  = y_output(:,4);

    model.fT3B    = y_output(:,5);
    model.T3TBGB  = y_output(:,6);
    model.T3TTRB  = y_output(:,7);
    model.T3ALBB  = y_output(:,8);

    model.TBGB    = y_output(:,9);
    model.TTRB    = y_output(:,10);
    model.ALBB    = y_output(:,11);

    model.fT4T    = y_output(:,12);
    model.T4TBGT  = y_output(:,13);
    model.T4TTRT  = y_output(:,14);
    model.T4ALBT  = y_output(:,15);
    model.fT3T    = y_output(:,16);
    model.T3TBGT  = y_output(:,17);
    model.T3TTRT  = y_output(:,18);
    model.T3ALBT  = y_output(:,19);
    model.TBGT    = y_output(:,20);
    model.TTRT    = y_output(:,21);
    model.ALBT    = y_output(:,22);

    model.fT4RB    = y_output(:,23);  
    model.T4TBGRB  = y_output(:,24);
    model.T4TTRRB  = y_output(:,25);
    model.T4ALBRB  = y_output(:,26);
    model.fT3RB    = y_output(:,27);
    model.T3TBGRB  = y_output(:,28);
    model.T3TTRRB  = y_output(:,29);
    model.T3ALBRB  = y_output(:,30);
    model.TBGRB    = y_output(:,31);
    model.TTRRB    = y_output(:,32);
    model.ALBRB    = y_output(:,33);
    model.T4RBT    = y_output(:,34);
    model.T3RBT    = y_output(:,35);

    model.fT4L    = y_output(:,36);
    model.T4TBGL  = y_output(:,37);
    model.T4TTRL  = y_output(:,38);
    model.T4ALBL  = y_output(:,39);
    model.fT3L    = y_output(:,40);
    model.T3TBGL  = y_output(:,41);
    model.T3TTRL  = y_output(:,42);
    model.T3ALBL  = y_output(:,43);
    model.T4LT    = y_output(:,44);
    model.T3LT    = y_output(:,45);

    model.XB      = y_output(:,46);
    model.XTTRB   = y_output(:,47);
    model.XT      = y_output(:,48);
    model.XTTRT   = y_output(:,49);
    model.XRB     = y_output(:,50);
    model.XTTRRB  = y_output(:,51);
    model.XL      = y_output(:,52);
    model.XTTRL   = y_output(:,53);
    model.XTBGB   = y_output(:,54);
    model.XTBGT   = y_output(:,55);
    model.XTBGRB  = y_output(:,56);
    model.XTBGL   = y_output(:,57);

    model.TBGL = (param.TBGtot - (model.TBGB + model.T4TBGB + model.T3TBGB + model.XTBGB)*param.VB - (model.TBGT + model.T4TBGT + model.T3TBGT + model.XTBGT)*param.VTB - (model.TBGRB + model.T4TBGRB + model.T3TBGRB + model.XTBGRB)*param.VRBB - (model.T4TBGL + model.T3TBGL + model.XTBGL)*param.VLB)/param.VLB;
    model.TTRL = (param.TTRtot - (model.TTRB + model.T4TTRB + model.T3TTRB + model.XTTRB)*param.VB - (model.TTRT + model.T4TTRT + model.T3TTRT + model.XTTRT)*param.VTB - (model.TTRRB + model.T4TTRRB + model.T3TTRRB + model.XTTRRB)*param.VRBB - (model.T4TTRL + model.T3TTRL + model.XTTRL)*param.VLB)/param.VLB;
    model.ALBL = (param.ALBtot - (model.ALBB + model.T4ALBB + model.T3ALBB)*param.VB - (model.ALBT + model.T4ALBT + model.T3ALBT)*param.VTB - (model.ALBRB + model.T4ALBRB + model.T3ALBRB)*param.VRBB - (model.T4ALBL + model.T3ALBL)*param.VLB)/param.VLB;



    figure(121)
    plot(hour, model.XB)
    hold on
    xlim([50*24 52*24])
    xlabel('Time (hour)')
    ylabel('Free X')

    figure(122)
    plot(hour, model.fT4B)
    hold on
    xlim([50*24 52*24])
    xlabel('Time (hour)')
    ylabel('fT4')

    figure(123)
    plot(hour, model.T4LT)
    hold on
    xlim([50*24 52*24])
    xlabel('Time (hour)')
    ylabel('Liver T4')

    figure(124)
    plot(hour, model.T4RBT)
    hold on
    xlim([50*24 52*24])
    xlabel('Time (hour)')
    ylabel('RB T4')

    figure(125)
    plot(hour, model.fT3B)
    hold on
    xlim([50*24 52*24])
    xlabel('Time (hour)')
    ylabel('fT3')

    figure(126)
    plot(hour, model.T3LT)
    hold on
    xlim([50*24 52*24])
    xlabel('Time (hour)')
    ylabel('Liver T3')

    figure(127)
    plot(hour, model.T3RBT)
    hold on
    xlim([50*24 52*24])
    xlabel('Time (hour)')
    ylabel('RB T3')

end


%%
toc
%

