function dydt = TH_PBK_Nonspatial_ODE_clamped(t, y, option, param)

%% -------------------------- PARAMETERS MAPPING ----------------------------------%%
k1 = param.k1;
k2 = param.k2;
k3 = param.k3;
k4 = param.k4;
k5 = param.k5;
k6 = param.k6;
k7 = param.k7;
k8 = param.k8;
k9 = param.k9;
k10 = param.k10;
k11 = param.k11;
k12 = param.k12;
k20 = param.k20;
k21 = param.k21;
k22 = param.k22;
k23 = param.k23;
k24 = param.k24;
k25 = param.k25;
k26 = param.k26;
k27 = param.k27;
k28 = param.k28;
k29 = param.k29;
k30 = param.k30;
k31 = param.k31;
k32 = param.k32;
k33 = param.k33;
k34 = param.k34;
k35 = param.k35;
k36 = param.k36;
k37 = param.k37;
k38 = param.k38;
k39 = param.k39;
k42 = param.k42;
k43 = param.k43;
fuT4RBT = param.fuT4RBT;
fuT3RBT = param.fuT3RBT;
fuT4LT = param.fuT4LT;
fuT3LT = param.fuT3LT;

TBGtot = param.TBGtot;
TTRtot = param.TTRtot;
ALBtot = param.ALBtot;

QC   = param.QC;
QT   = param.QT;
QRB  = param.QRB;
QL   = param.QL;
VB   = param.VB;
VTB  = param.VTB;
VRBB = param.VRBB;
VRBT = param.VRBT;
VLB  = param.VLB;
VLT  = param.VLT;


%% ------------------------- STATE NAME MAPPING----------------------------%%
fT4B    = y(1);
T4TBGB  = y(2);
T4TTRB  = y(3);
T4ALBB  = y(4);
fT3B    = y(5);
T3TBGB  = y(6);
T3TTRB  = y(7);
T3ALBB  = y(8);
TBGB    = y(9);
TTRB    = y(10);
ALBB    = y(11);

fT4T    = y(12);
T4TBGT  = y(13);
T4TTRT  = y(14);
T4ALBT  = y(15);
fT3T    = y(16);
T3TBGT  = y(17);
T3TTRT  = y(18);
T3ALBT  = y(19);
TBGT    = y(20);
TTRT    = y(21);
ALBT    = y(22);

fT4RB    = y(23);
T4TBGRB  = y(24);
T4TTRRB  = y(25);
T4ALBRB  = y(26);
fT3RB    = y(27);
T3TBGRB  = y(28);
T3TTRRB  = y(29);
T3ALBRB  = y(30);
TBGRB    = y(31);
TTRRB    = y(32);
ALBRB    = y(33);
T4RBT    = y(34);
T3RBT    = y(35);

fT4L    = y(36);
T4TBGL  = y(37);
T4TTRL  = y(38);
T4ALBL  = y(39);
fT3L    = y(40);
T3TBGL  = y(41);
T3TTRL  = y(42);
T3ALBL  = y(43);
T4LT    = y(44);
T3LT    = y(45);

XB      = y(46);
XTTRB   = y(47);
XT      = y(48);
XTTRT   = y(49);
XRB     = y(50);
XTTRRB  = y(51);
XL      = y(52);
XTTRL   = y(53);
XTBGB   = y(54);
XTBGT   = y(55);
XTBGRB  = y(56);
XTBGL   = y(57);

TBGL = (TBGtot - (TBGB + T4TBGB + T3TBGB + XTBGB)*VB - (TBGT + T4TBGT + T3TBGT + XTBGT)*VTB - (TBGRB + T4TBGRB + T3TBGRB + XTBGRB)*VRBB - (T4TBGL + T3TBGL + XTBGL)*VLB)/VLB;
TTRL = (TTRtot - (TTRB + T4TTRB + T3TTRB + XTTRB)*VB - (TTRT + T4TTRT + T3TTRT + XTTRT)*VTB - (TTRRB + T4TTRRB + T3TTRRB + XTTRRB)*VRBB - (T4TTRL + T3TTRL + XTTRL)*VLB)/VLB;
ALBL = (ALBtot - (ALBB + T4ALBB + T3ALBB)*VB - (ALBT + T4ALBT + T3ALBT)*VTB - (ALBRB + T4ALBRB + T3ALBRB)*VRBB - (T4ALBL + T3ALBL)*VLB)/VLB;


%% ------------------------------ ODEs-------------------------------------%%

dydt = zeros(length(y),1); %make dydt as a column vector as required by MatLab ode function

%Body Blood Compartment 
%fT4B 
dydt(1) = 0;%-k1*fT4B*TBGB + k2*T4TBGB - k3*fT4B*TTRB + k4*T4TTRB - k5*fT4B*ALBB + k6*T4ALBB + (fT4T*QT + fT4RB*QRB + fT4L*QL - fT4B*QC)/VB;

%T4TBGB
dydt(2)	= k1*fT4B*TBGB - k2*T4TBGB + (T4TBGT*QT + T4TBGRB*QRB + T4TBGL*QL - T4TBGB*QC)/VB;

%T4TTRB
dydt(3) = k3*fT4B*TTRB - k4*T4TTRB + (T4TTRT*QT + T4TTRRB*QRB + T4TTRL*QL - T4TTRB*QC)/VB;

%T4ALBB
dydt(4) = k5*fT4B*ALBB - k6*T4ALBB + (T4ALBT*QT + T4ALBRB*QRB + T4ALBL*QL - T4ALBB*QC)/VB;

%fT3B
dydt(5) = 0;%-k7*fT3B*TBGB + k8*T3TBGB - k9*fT3B*TTRB + k10*T3TTRB - k11*fT3B*ALBB + k12*T3ALBB + (fT3T*QT + fT3RB*QRB + fT3L*QL - fT3B*QC)/VB;

%T3TBGB
dydt(6)	= k7*fT3B*TBGB - k8*T3TBGB + (T3TBGT*QT + T3TBGRB*QRB + T3TBGL*QL - T3TBGB*QC)/VB;

%T3TTRB
dydt(7) = k9*fT3B*TTRB - k10*T3TTRB + (T3TTRT*QT + T3TTRRB*QRB + T3TTRL*QL - T3TTRB*QC)/VB;

%T3ALBB
dydt(8) = k11*fT3B*ALBB - k12*T3ALBB + (T3ALBT*QT + T3ALBRB*QRB + T3ALBL*QL - T3ALBB*QC)/VB;

%TBGB
dydt(9) = -k1*fT4B*TBGB + k2*T4TBGB - k7*fT3B*TBGB + k8*T3TBGB - k42*XB*TBGB + k43*XTBGB + (TBGT*QT + TBGRB*QRB + TBGL*QL - TBGB*QC)/VB;

%TTRB
dydt(10) = -k3*fT4B*TTRB + k4*T4TTRB - k9*fT3B*TTRB + k10*T3TTRB - k38*XB*TTRB + k39*XTTRB + (TTRT*QT + TTRRB*QRB + TTRL*QL - TTRB*QC)/VB;

%ALBB
dydt(11) = -k5*fT4B*ALBB + k6*T4ALBB - k11*fT3B*ALBB + k12*T3ALBB + (ALBT*QT + ALBRB*QRB + ALBL*QL - ALBB*QC)/VB;


%Thyroid blood Compartment
%fT4T
dydt(12) = -k1*fT4T*TBGT + k2*T4TBGT - k3*fT4T*TTRT + k4*T4TTRT - k5*fT4T*ALBT + k6*T4ALBT + k20/VTB + (fT4B-fT4T)*QT/VTB;

%T4TBGT
dydt(13) = k1*fT4T*TBGT - k2*T4TBGT + (T4TBGB-T4TBGT)*QT/VTB;

%T4TTRT
dydt(14) = k3*fT4T*TTRT - k4*T4TTRT + (T4TTRB-T4TTRT)*QT/VTB;

%T4ALBT
dydt(15) = k5*fT4T*ALBT - k6*T4ALBT + (T4ALBB-T4ALBT)*QT/VTB;

%fT3T
dydt(16) = -k7*fT3T*TBGT + k8*T3TBGT - k9*fT3T*TTRT + k10*T3TTRT - k11*fT3T*ALBT + k12*T3ALBT + k22/VTB + (fT3B-fT3T)*QT/VTB;

%T3TBGT
dydt(17) = k7*fT3T*TBGT - k8*T3TBGT + (T3TBGB-T3TBGT)*QT/VTB;

%T3TTRT
dydt(18) = k9*fT3T*TTRT - k10*T3TTRT + (T3TTRB-T3TTRT)*QT/VTB;

%T3ALBT
dydt(19) = k11*fT3T*ALBT - k12*T3ALBT + (T3ALBB-T3ALBT)*QT/VTB;

%TBGT
dydt(20) = -k1*fT4T*TBGT + k2*T4TBGT - k7*fT3T*TBGT + k8*T3TBGT - k42*XT*TBGT + k43*XTBGT + (TBGB-TBGT)*QT/VTB;

%TTRT
dydt(21) = -k3*fT4T*TTRT + k4*T4TTRT - k9*fT3T*TTRT + k10*T3TTRT - k38*XT*TTRT + k39*XTTRT + (TTRB-TTRT)*QT/VTB;

%ALBT
dydt(22) = -k5*fT4T*ALBT + k6*T4ALBT - k11*fT3T*ALBT + k12*T3ALBT + (ALBB-ALBT)*QT/VTB;


%RB Compartment
%fT4RB
dydt(23) = -k1*fT4RB*TBGRB + k2*T4TBGRB - k3*fT4RB*TTRRB + k4*T4TTRRB - k5*fT4RB*ALBRB + k6*T4ALBRB + (-k21*fT4RB + k28*T4RBT*fuT4RBT)/VRBB + (fT4B-fT4RB)*QRB/VRBB;

%T4TBGRB
dydt(24) = k1*fT4RB*TBGRB - k2*T4TBGRB + (T4TBGB-T4TBGRB)*QRB/VRBB;

%T4TTRRB
dydt(25) = k3*fT4RB*TTRRB - k4*T4TTRRB + (T4TTRB-T4TTRRB)*QRB/VRBB;

%T4ALBRB
dydt(26) = k5*fT4RB*ALBRB - k6*T4ALBRB + (T4ALBB-T4ALBRB)*QRB/VRBB;

%fT3RB
dydt(27) = -k7*fT3RB*TBGRB + k8*T3TBGRB - k9*fT3RB*TTRRB + k10*T3TTRRB - k11*fT3RB*ALBRB + k12*T3ALBRB + (-k23*fT3RB + k29*T3RBT*fuT3RBT)/VRBB + (fT3B-fT3RB)*QRB/VRBB;

%T3TBGRB
dydt(28) = k7*fT3RB*TBGRB - k8*T3TBGRB + (T3TBGB-T3TBGRB)*QRB/VRBB;

%T3TTRRB
dydt(29) = k9*fT3RB*TTRRB - k10*T3TTRRB + (T3TTRB-T3TTRRB)*QRB/VRBB;

%T3ALBRB
dydt(30) = k11*fT3RB*ALBRB - k12*T3ALBRB + (T3ALBB-T3ALBRB)*QRB/VRBB;

%TBGRB
dydt(31) = -k1*fT4RB*TBGRB + k2*T4TBGRB - k7*fT3RB*TBGRB + k8*T3TBGRB - k42*XRB*TBGRB + k43*XTBGRB + (TBGB-TBGRB)*QRB/VRBB;

%TTRRB
dydt(32) = -k3*fT4RB*TTRRB + k4*T4TTRRB - k9*fT3RB*TTRRB + k10*T3TTRRB - k38*XRB*TTRRB + k39*XTTRRB + (TTRB-TTRRB)*QRB/VRBB;

%ALBRB
dydt(33) = -k5*fT4RB*ALBRB + k6*T4ALBRB - k11*fT3RB*ALBRB + k12*T3ALBRB + (ALBB-ALBRB)*QRB/VRBB;

%T4RBT
dydt(34) = (k21*fT4RB - k28*T4RBT*fuT4RBT)/VRBT - k24*T4RBT*fuT4RBT - k32*T4RBT*fuT4RBT;

%T3RBT
dydt(35) = (k23*fT3RB - k29*T3RBT*fuT3RBT)/VRBT + k24*T4RBT*fuT4RBT - k33*T3RBT*fuT3RBT;


%Liver Compartment
%fT4L
dydt(36) = -k1*fT4L*TBGL + k2*T4TBGL - k3*fT4L*TTRL + k4*T4TTRL - k5*fT4L*ALBL + k6*T4ALBL + (-k25*fT4L + k30*T4LT*fuT4LT)/VLB + (fT4B-fT4L)*QL/VLB;

%T4TBGL
dydt(37) = k1*fT4L*TBGL - k2*T4TBGL + (T4TBGB-T4TBGL)*QL/VLB;

%T4TTRL
dydt(38) = k3*fT4L*TTRL - k4*T4TTRL + (T4TTRB-T4TTRL)*QL/VLB;

%T4ALBL
dydt(39) = k5*fT4L*ALBL - k6*T4ALBL + (T4ALBB-T4ALBL)*QL/VLB;

%fT3L
dydt(40) = -k7*fT3L*TBGL + k8*T3TBGL - k9*fT3L*TTRL + k10*T3TTRL - k11*fT3L*ALBL + k12*T3ALBL + (-k27*fT3L + k31*T3LT*fuT3LT)/VLB + (fT3B-fT3L)*QL/VLB;

%T3TBGL
dydt(41) = k7*fT3L*TBGL - k8*T3TBGL + (T3TBGB-T3TBGL)*QL/VLB;

%T3TTRL
dydt(42) = k9*fT3L*TTRL - k10*T3TTRL + (T3TTRB-T3TTRL)*QL/VLB;

%T3ALBL
dydt(43) = k11*fT3L*ALBL - k12*T3ALBL + (T3ALBB-T3ALBL)*QL/VLB;

%T4LT
dydt(44) = (k25*fT4L - k30*T4LT*fuT4LT)/VLT - k26*T4LT*fuT4LT - k34*T4LT*fuT4LT;

%T3LT
dydt(45) = (k27*fT3L - k31*T3LT*fuT3LT)/VLT + k26*T4LT*fuT4LT - k35*T3LT*fuT3LT;


%Binding Between X and TTR
%XB (Body Blood Compartment)
dydt(46) = k36/VB - k37*XB - k38*XB*TTRB - k42*XB*TBGB + k39*XTTRB + k43*XTBGB + (XT*QT + XRB*QRB + XL*QL - XB*QC)/VB;

%XTTRB (Body Blood Compartment)
dydt(47) = k38*XB*TTRB - k39*XTTRB + (XTTRT*QT + XTTRRB*QRB + XTTRL*QL - XTTRB*QC)/VB;

%XT (Thyroid blood Compartment)
dydt(48) = -k38*XT*TTRT - k42*XT*TBGT + k39*XTTRT + k43*XTBGT + (XB-XT)*QT/VTB;

%XTTRT (Thyroid blood Compartment)
dydt(49) = k38*XT*TTRT - k39*XTTRT + (XTTRB-XTTRT)*QT/VTB;

%XRB (RB blood Compartment)
dydt(50) = -k38*XRB*TTRRB - k42*XRB*TBGRB + k39*XTTRRB + k43*XTBGRB + (XB-XRB)*QRB/VRBB;

%XTTRRB (RB blood Compartment)
dydt(51) = k38*XRB*TTRRB - k39*XTTRRB + (XTTRB-XTTRRB)*QRB/VRBB;

%XL (Liver Blood Compartment)
dydt(52) = -k38*XL*TTRL - k42*XL*TBGL + k39*XTTRL + k43*XTBGL + (XB-XL)*QL/VLB;

%XTTRL (Liver blood Compartment)
dydt(53) = k38*XL*TTRL - k39*XTTRL + (XTTRB-XTTRL)*QL/VLB;

%Binding Between X and TBG
%XTBGB (Body Blood Compartment)
dydt(54) = k42*XB*TBGB - k43*XTBGB + (XTBGT*QT + XTBGRB*QRB + XTBGL*QL - XTBGB*QC)/VB;

%XTBGT (Thyroid blood Compartment)
dydt(55) = k42*XT*TBGT - k43*XTBGT + (XTBGB-XTBGT)*QT/VTB;

%XTBGRB (RB blood Compartment)
dydt(56) = k42*XRB*TBGRB - k43*XTBGRB + (XTBGB-XTBGRB)*QRB/VRBB;

%XTBGL (Liver blood Compartment)
dydt(57) = k42*XL*TBGL - k43*XTBGL + (XTBGB-XTBGL)*QL/VLB;

end