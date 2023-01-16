function dydt = TH_PBK_Nonspatial_ODE_tracer(t, y, option, param)

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
tfT4B   = y(1);
tT4TBGB  = y(2);
tT4TTRB  = y(3);
tT4ALBB  = y(4);
fT4B    = y(5);
T4TBGB  = y(6);
T4TTRB  = y(7);
T4ALBB  = y(8);
tfT3B    = y(9);
tT3TBGB  = y(10);
tT3TTRB  = y(11);
tT3ALBB  = y(12);
fT3B    = y(13);
T3TBGB  = y(14);
T3TTRB  = y(15);
T3ALBB  = y(16);
TBGB    = y(17);
TTRB    = y(18);
ALBB    = y(19);

tfT4T   = y(20);
tT4TBGT  = y(21);
tT4TTRT  = y(22);
tT4ALBT  = y(23);
fT4T    = y(24);
T4TBGT  = y(25);
T4TTRT  = y(26);
T4ALBT  = y(27);
tfT3T    = y(28);
tT3TBGT  = y(29);
tT3TTRT  = y(30);
tT3ALBT  = y(31);
fT3T    = y(32);
T3TBGT  = y(33);
T3TTRT  = y(34);
T3ALBT  = y(35);
TBGT    = y(36);
TTRT    = y(37);
ALBT    = y(38);

tfT4RB   = y(39);
tT4TBGRB  = y(40);
tT4TTRRB  = y(41);
tT4ALBRB  = y(42);
fT4RB    = y(43);
T4TBGRB  = y(44);
T4TTRRB  = y(45);
T4ALBRB  = y(46);
tfT3RB    = y(47);
tT3TBGRB  = y(48);
tT3TTRRB  = y(49);
tT3ALBRB  = y(50);
fT3RB    = y(51);
T3TBGRB  = y(52);
T3TTRRB  = y(53);
T3ALBRB  = y(54);
TBGRB    = y(55);
TTRRB    = y(56);
ALBRB    = y(57);
tT4RBT    = y(58);
tT3RBT    = y(59);
T4RBT    = y(60);
T3RBT    = y(61);

tfT4L   = y(62);
tT4TBGL  = y(63);
tT4TTRL  = y(64);
tT4ALBL  = y(65);
fT4L    = y(66);
T4TBGL  = y(67);
T4TTRL  = y(68);
T4ALBL  = y(69);
tfT3L    = y(70);
tT3TBGL  = y(71);
tT3TTRL  = y(72);
tT3ALBL  = y(73);
fT3L    = y(74);
T3TBGL  = y(75);
T3TTRL  = y(76);
T3ALBL  = y(77);
tT4LT    = y(78);
tT3LT    = y(79);
T4LT    = y(80);
T3LT    = y(81);

XB      = y(82);
XTTRB   = y(83);
XT      = y(84);
XTTRT   = y(85);
XRB     = y(86);
XTTRRB  = y(87);
XL      = y(88);
XTTRL   = y(89);
XTBGB   = y(90);
XTBGT   = y(91);
XTBGRB  = y(92);
XTBGL   = y(93);

TBGL = (TBGtot - (TBGB + T4TBGB + T3TBGB + tT4TBGB + tT3TBGB + XTBGB)*VB - (TBGT + T4TBGT + T3TBGT + tT4TBGT + tT3TBGT + XTBGT)*VTB - (TBGRB + T4TBGRB + T3TBGRB + tT4TBGRB + tT3TBGRB + XTBGRB)*VRBB - (T4TBGL + T3TBGL + tT4TBGL + tT3TBGL + XTBGL)*VLB)/VLB;
TTRL = (TTRtot - (TTRB + T4TTRB + T3TTRB + tT4TTRB + tT3TTRB + XTTRB)*VB - (TTRT + T4TTRT + T3TTRT + tT4TTRT + tT3TTRT + XTTRT)*VTB - (TTRRB + T4TTRRB + T3TTRRB + tT4TTRRB + tT3TTRRB + XTTRRB)*VRBB - (T4TTRL + T3TTRL + tT4TTRL + tT3TTRL + XTTRL)*VLB)/VLB;
ALBL = (ALBtot - (ALBB + T4ALBB + T3ALBB + tT4ALBB + tT3ALBB)*VB - (ALBT + T4ALBT + T3ALBT + tT4ALBT + tT3ALBT)*VTB - (ALBRB + T4ALBRB + T3ALBRB + tT4ALBRB + tT3ALBRB)*VRBB - (T4ALBL + T3ALBL + tT4ALBL + tT3ALBL)*VLB)/VLB;


%% ------------------------------ ODEs-------------------------------------%%
dydt = zeros(length(y),1); %make dydt as a column vector as required by MatLab ode function

%Body Blood Compartment 
%tfT4B (free T4 tracer)
dydt(1) = -k1*tfT4B*TBGB + k2*tT4TBGB - k3*tfT4B*TTRB + k4*tT4TTRB - k5*tfT4B*ALBB + k6*tT4ALBB + (tfT4T*QT + tfT4RB*QRB + tfT4L*QL - tfT4B*QC)/VB;

%tT4TBGB (T4 tracer TBG)
dydt(2)	= k1*tfT4B*TBGB - k2*tT4TBGB + (tT4TBGT*QT + tT4TBGRB*QRB + tT4TBGL*QL - tT4TBGB*QC)/VB;

%tT4TTRB (T4 tracer TTR)
dydt(3) = k3*tfT4B*TTRB - k4*tT4TTRB + (tT4TTRT*QT + tT4TTRRB*QRB + tT4TTRL*QL - tT4TTRB*QC)/VB;

%tT4ALBB (T4 tracer ALB)
dydt(4) = k5*tfT4B*ALBB - k6*tT4ALBB + (tT4ALBT*QT + tT4ALBRB*QRB + tT4ALBL*QL - tT4ALBB*QC)/VB;

%fT4B 
dydt(5) = -k1*fT4B*TBGB + k2*T4TBGB - k3*fT4B*TTRB + k4*T4TTRB - k5*fT4B*ALBB + k6*T4ALBB + (fT4T*QT + fT4RB*QRB + fT4L*QL - fT4B*QC)/VB;

%T4TBGB
dydt(6)	= k1*fT4B*TBGB - k2*T4TBGB + (T4TBGT*QT + T4TBGRB*QRB + T4TBGL*QL - T4TBGB*QC)/VB;

%T4TTRB
dydt(7) = k3*fT4B*TTRB - k4*T4TTRB + (T4TTRT*QT + T4TTRRB*QRB + T4TTRL*QL - T4TTRB*QC)/VB;

%T4ALBB
dydt(8) = k5*fT4B*ALBB - k6*T4ALBB + (T4ALBT*QT + T4ALBRB*QRB + T4ALBL*QL - T4ALBB*QC)/VB;

%tfT3B (free T3 tracer)
dydt(9) = -k7*tfT3B*TBGB + k8*tT3TBGB - k9*tfT3B*TTRB + k10*tT3TTRB - k11*tfT3B*ALBB + k12*tT3ALBB + (tfT3T*QT + tfT3RB*QRB + tfT3L*QL - tfT3B*QC)/VB;

%tT3TBGB (T3 tracer TBG)
dydt(10) = k7*tfT3B*TBGB - k8*tT3TBGB + (tT3TBGT*QT + tT3TBGRB*QRB + tT3TBGL*QL - tT3TBGB*QC)/VB;

%tT3TTRB (T3 tracer TTR)
dydt(11) = k9*tfT3B*TTRB - k10*tT3TTRB + (tT3TTRT*QT + tT3TTRRB*QRB + tT3TTRL*QL - tT3TTRB*QC)/VB;

%tT3ALBB (T3 tracer ALB)
dydt(12) = k11*tfT3B*ALBB - k12*tT3ALBB + (tT3ALBT*QT + tT3ALBRB*QRB + tT3ALBL*QL - tT3ALBB*QC)/VB;

%fT3B
dydt(13) = -k7*fT3B*TBGB + k8*T3TBGB - k9*fT3B*TTRB + k10*T3TTRB - k11*fT3B*ALBB + k12*T3ALBB + (fT3T*QT + fT3RB*QRB + fT3L*QL - fT3B*QC)/VB;

%T3TBGB
dydt(14) = k7*fT3B*TBGB - k8*T3TBGB + (T3TBGT*QT + T3TBGRB*QRB + T3TBGL*QL - T3TBGB*QC)/VB;

%T3TTRB
dydt(15) = k9*fT3B*TTRB - k10*T3TTRB + (T3TTRT*QT + T3TTRRB*QRB + T3TTRL*QL - T3TTRB*QC)/VB;

%T3ALBB
dydt(16) = k11*fT3B*ALBB - k12*T3ALBB + (T3ALBT*QT + T3ALBRB*QRB + T3ALBL*QL - T3ALBB*QC)/VB;

%TBGB
dydt(17) = -k1*fT4B*TBGB + k2*T4TBGB - k1*tfT4B*TBGB + k2*tT4TBGB - k7*fT3B*TBGB + k8*T3TBGB - k7*tfT3B*TBGB + k8*tT3TBGB - k42*XB*TBGB + k43*XTBGB + (TBGT*QT + TBGRB*QRB + TBGL*QL - TBGB*QC)/VB;

%TTRB
dydt(18) = -k3*fT4B*TTRB + k4*T4TTRB - k3*tfT4B*TTRB + k4*tT4TTRB - k9*fT3B*TTRB + k10*T3TTRB - k9*tfT3B*TTRB + k10*tT3TTRB - k38*XB*TTRB + k39*XTTRB + (TTRT*QT + TTRRB*QRB + TTRL*QL - TTRB*QC)/VB;

%ALBB
dydt(19) = -k5*fT4B*ALBB + k6*T4ALBB - k5*tfT4B*ALBB + k6*tT4ALBB - k11*fT3B*ALBB + k12*T3ALBB - k11*tfT3B*ALBB + k12*tT3ALBB + (ALBT*QT + ALBRB*QRB + ALBL*QL - ALBB*QC)/VB;


%Thyroid blood Compartment
%tfT4T
dydt(20) = -k1*tfT4T*TBGT + k2*tT4TBGT - k3*tfT4T*TTRT + k4*tT4TTRT - k5*tfT4T*ALBT + k6*tT4ALBT + (tfT4B-tfT4T)*QT/VTB;

%tT4TBGT
dydt(21) = k1*tfT4T*TBGT - k2*tT4TBGT + (tT4TBGB-tT4TBGT)*QT/VTB;

%tT4TTRT
dydt(22) = k3*tfT4T*TTRT - k4*tT4TTRT + (tT4TTRB-tT4TTRT)*QT/VTB;

%tT4ALBT
dydt(23) = k5*tfT4T*ALBT - k6*tT4ALBT + (tT4ALBB-tT4ALBT)*QT/VTB;

%fT4T
dydt(24) = -k1*fT4T*TBGT + k2*T4TBGT - k3*fT4T*TTRT + k4*T4TTRT - k5*fT4T*ALBT + k6*T4ALBT + k20/VTB + (fT4B-fT4T)*QT/VTB;

%T4TBGT
dydt(25) = k1*fT4T*TBGT - k2*T4TBGT + (T4TBGB-T4TBGT)*QT/VTB;

%T4TTRT
dydt(26) = k3*fT4T*TTRT - k4*T4TTRT + (T4TTRB-T4TTRT)*QT/VTB;

%T4ALBT
dydt(27) = k5*fT4T*ALBT - k6*T4ALBT + (T4ALBB-T4ALBT)*QT/VTB;

%tfT3T
dydt(28) = -k7*tfT3T*TBGT + k8*tT3TBGT - k9*tfT3T*TTRT + k10*tT3TTRT - k11*tfT3T*ALBT + k12*tT3ALBT + (tfT3B-tfT3T)*QT/VTB;

%tT3TBGT
dydt(29) = k7*tfT3T*TBGT - k8*tT3TBGT + (tT3TBGB-tT3TBGT)*QT/VTB;

%tT3TTRT
dydt(30) = k9*tfT3T*TTRT - k10*tT3TTRT + (tT3TTRB-tT3TTRT)*QT/VTB;

%tT3ALBT
dydt(31) = k11*tfT3T*ALBT - k12*tT3ALBT + (tT3ALBB-tT3ALBT)*QT/VTB;

%fT3T
dydt(32) = -k7*fT3T*TBGT + k8*T3TBGT - k9*fT3T*TTRT + k10*T3TTRT - k11*fT3T*ALBT + k12*T3ALBT + k22/VTB + (fT3B-fT3T)*QT/VTB;

%T3TBGT
dydt(33) = k7*fT3T*TBGT - k8*T3TBGT + (T3TBGB-T3TBGT)*QT/VTB;

%T3TTRT
dydt(34) = k9*fT3T*TTRT - k10*T3TTRT + (T3TTRB-T3TTRT)*QT/VTB;

%T3ALBT
dydt(35) = k11*fT3T*ALBT - k12*T3ALBT + (T3ALBB-T3ALBT)*QT/VTB;

%TBGT
dydt(36) = -k1*fT4T*TBGT + k2*T4TBGT - k1*tfT4T*TBGT + k2*tT4TBGT - k7*fT3T*TBGT + k8*T3TBGT - k7*tfT3T*TBGT + k8*tT3TBGT - k42*XT*TBGT + k43*XTBGT + (TBGB-TBGT)*QT/VTB;

%TTRT
dydt(37) = -k3*fT4T*TTRT + k4*T4TTRT - k3*tfT4T*TTRT + k4*tT4TTRT - k9*fT3T*TTRT + k10*T3TTRT - k9*tfT3T*TTRT + k10*tT3TTRT - k38*XT*TTRT + k39*XTTRT + (TTRB-TTRT)*QT/VTB;

%ALBT
dydt(38) = -k5*fT4T*ALBT + k6*T4ALBT - k5*tfT4T*ALBT + k6*tT4ALBT - k11*fT3T*ALBT + k12*T3ALBT - k11*tfT3T*ALBT + k12*tT3ALBT + (ALBB-ALBT)*QT/VTB;


%RB Compartment
%tfT4RB
dydt(39) = -k1*tfT4RB*TBGRB + k2*tT4TBGRB - k3*tfT4RB*TTRRB + k4*tT4TTRRB - k5*tfT4RB*ALBRB + k6*tT4ALBRB + (-k21*tfT4RB + k28*tT4RBT*fuT4RBT)/VRBB + (tfT4B-tfT4RB)*QRB/VRBB;

%tT4TBGRB
dydt(40) = k1*tfT4RB*TBGRB - k2*tT4TBGRB + (tT4TBGB-tT4TBGRB)*QRB/VRBB;

%tT4TTRRB
dydt(41) = k3*tfT4RB*TTRRB - k4*tT4TTRRB + (tT4TTRB-tT4TTRRB)*QRB/VRBB;

%tT4ALBRB
dydt(42) = k5*tfT4RB*ALBRB - k6*tT4ALBRB + (tT4ALBB-tT4ALBRB)*QRB/VRBB;

%fT4RB
dydt(43) = -k1*fT4RB*TBGRB + k2*T4TBGRB - k3*fT4RB*TTRRB + k4*T4TTRRB - k5*fT4RB*ALBRB + k6*T4ALBRB + (-k21*fT4RB + k28*T4RBT*fuT4RBT)/VRBB + (fT4B-fT4RB)*QRB/VRBB;

%T4TBGRB
dydt(44) = k1*fT4RB*TBGRB - k2*T4TBGRB + (T4TBGB-T4TBGRB)*QRB/VRBB;

%T4TTRRB
dydt(45) = k3*fT4RB*TTRRB - k4*T4TTRRB + (T4TTRB-T4TTRRB)*QRB/VRBB;

%T4ALBRB
dydt(46) = k5*fT4RB*ALBRB - k6*T4ALBRB + (T4ALBB-T4ALBRB)*QRB/VRBB;

%tfT3RB
dydt(47) = -k7*tfT3RB*TBGRB + k8*tT3TBGRB - k9*tfT3RB*TTRRB + k10*tT3TTRRB - k11*tfT3RB*ALBRB + k12*tT3ALBRB + (-k23*tfT3RB + k29*tT3RBT*fuT3RBT)/VRBB + (tfT3B-tfT3RB)*QRB/VRBB;

%tT3TBGRB
dydt(48) = k7*tfT3RB*TBGRB - k8*tT3TBGRB + (tT3TBGB-tT3TBGRB)*QRB/VRBB;

%tT3TTRRB
dydt(49) = k9*tfT3RB*TTRRB - k10*tT3TTRRB + (tT3TTRB-tT3TTRRB)*QRB/VRBB;

%tT3ALBRB
dydt(50) = k11*tfT3RB*ALBRB - k12*tT3ALBRB + (tT3ALBB-tT3ALBRB)*QRB/VRBB;

%fT3RB
dydt(51) = -k7*fT3RB*TBGRB + k8*T3TBGRB - k9*fT3RB*TTRRB + k10*T3TTRRB - k11*fT3RB*ALBRB + k12*T3ALBRB + (-k23*fT3RB + k29*T3RBT*fuT3RBT)/VRBB + (fT3B-fT3RB)*QRB/VRBB;

%T3TBGRB
dydt(52) = k7*fT3RB*TBGRB - k8*T3TBGRB + (T3TBGB-T3TBGRB)*QRB/VRBB;

%T3TTRRB
dydt(53) = k9*fT3RB*TTRRB - k10*T3TTRRB + (T3TTRB-T3TTRRB)*QRB/VRBB;

%T3ALBRB
dydt(54) = k11*fT3RB*ALBRB - k12*T3ALBRB + (T3ALBB-T3ALBRB)*QRB/VRBB;

%TBGRB
dydt(55) = -k1*fT4RB*TBGRB + k2*T4TBGRB - k1*tfT4RB*TBGRB + k2*tT4TBGRB - k7*fT3RB*TBGRB + k8*T3TBGRB - k7*tfT3RB*TBGRB + k8*tT3TBGRB - k42*XRB*TBGRB + k43*XTBGRB + (TBGB-TBGRB)*QRB/VRBB;

%TTRRB
dydt(56) = -k3*fT4RB*TTRRB + k4*T4TTRRB - k3*tfT4RB*TTRRB + k4*tT4TTRRB - k9*fT3RB*TTRRB + k10*T3TTRRB - k9*tfT3RB*TTRRB + k10*tT3TTRRB - k38*XRB*TTRRB + k39*XTTRRB + (TTRB-TTRRB)*QRB/VRBB;

%ALBRB
dydt(57) = -k5*fT4RB*ALBRB + k6*T4ALBRB - k5*tfT4RB*ALBRB + k6*tT4ALBRB - k11*fT3RB*ALBRB + k12*T3ALBRB - k11*tfT3RB*ALBRB + k12*tT3ALBRB + (ALBB-ALBRB)*QRB/VRBB;

%tT4RBT
dydt(58) = (k21*tfT4RB - k28*tT4RBT*fuT4RBT)/VRBT - k24*tT4RBT*fuT4RBT - k32*tT4RBT*fuT4RBT;

%tT3RBT
dydt(59) = (k23*tfT3RB - k29*tT3RBT*fuT3RBT)/VRBT + k24*tT4RBT*fuT4RBT - k33*tT3RBT*fuT3RBT;

%T4RBT
dydt(60) = (k21*fT4RB - k28*T4RBT*fuT4RBT)/VRBT - k24*T4RBT*fuT4RBT - k32*T4RBT*fuT4RBT;

%T3RBT
dydt(61) = (k23*fT3RB - k29*T3RBT*fuT3RBT)/VRBT + k24*T4RBT*fuT4RBT  - k33*T3RBT*fuT3RBT;


%Liver Compartment
%tfT4L
dydt(62) = -k1*tfT4L*TBGL + k2*tT4TBGL - k3*tfT4L*TTRL + k4*tT4TTRL - k5*tfT4L*ALBL + k6*tT4ALBL + (-k25*tfT4L + k30*tT4LT*fuT4LT)/VLB + (tfT4B-tfT4L)*QL/VLB;

%tT4TBGL
dydt(63) = k1*tfT4L*TBGL - k2*tT4TBGL + (tT4TBGB-tT4TBGL)*QL/VLB;

%tT4TTRL
dydt(64) = k3*tfT4L*TTRL - k4*tT4TTRL + (tT4TTRB-tT4TTRL)*QL/VLB;

%tT4ALBL
dydt(65) = k5*tfT4L*ALBL - k6*tT4ALBL + (tT4ALBB-tT4ALBL)*QL/VLB;

%fT4L
dydt(66) = -k1*fT4L*TBGL + k2*T4TBGL - k3*fT4L*TTRL + k4*T4TTRL - k5*fT4L*ALBL + k6*T4ALBL + (-k25*fT4L + k30*T4LT*fuT4LT)/VLB + (fT4B-fT4L)*QL/VLB;

%T4TBGL
dydt(67) = k1*fT4L*TBGL - k2*T4TBGL + (T4TBGB-T4TBGL)*QL/VLB;

%T4TTRL
dydt(68) = k3*fT4L*TTRL - k4*T4TTRL + (T4TTRB-T4TTRL)*QL/VLB;

%T4ALBL
dydt(69) = k5*fT4L*ALBL - k6*T4ALBL + (T4ALBB-T4ALBL)*QL/VLB;

%tfT3L
dydt(70) = -k7*tfT3L*TBGL + k8*tT3TBGL - k9*tfT3L*TTRL + k10*tT3TTRL - k11*tfT3L*ALBL + k12*tT3ALBL + (-k27*tfT3L + k31*tT3LT*fuT3LT)/VLB + (tfT3B-tfT3L)*QL/VLB;

%tT3TBGL
dydt(71) = k7*tfT3L*TBGL - k8*tT3TBGL + (tT3TBGB-tT3TBGL)*QL/VLB;

%tT3TTRL
dydt(72) = k9*tfT3L*TTRL - k10*tT3TTRL + (tT3TTRB-tT3TTRL)*QL/VLB;

%tT3ALBL
dydt(73) = k11*tfT3L*ALBL - k12*tT3ALBL + (tT3ALBB-tT3ALBL)*QL/VLB;

%fT3L
dydt(74) = -k7*fT3L*TBGL + k8*T3TBGL - k9*fT3L*TTRL + k10*T3TTRL - k11*fT3L*ALBL + k12*T3ALBL + (-k27*fT3L + k31*T3LT*fuT3LT)/VLB + (fT3B-fT3L)*QL/VLB;

%T3TBGL
dydt(75) = k7*fT3L*TBGL - k8*T3TBGL + (T3TBGB-T3TBGL)*QL/VLB;

%T3TTRL
dydt(76) = k9*fT3L*TTRL - k10*T3TTRL + (T3TTRB-T3TTRL)*QL/VLB;

%T3ALBL
dydt(77) = k11*fT3L*ALBL - k12*T3ALBL + (T3ALBB-T3ALBL)*QL/VLB;

%tT4LT
dydt(78) = (k25*tfT4L - k30*tT4LT*fuT4LT)/VLT - k26*tT4LT*fuT4LT - k34*tT4LT*fuT4LT;

%tT3LT
dydt(79) = (k27*tfT3L - k31*tT3LT*fuT3LT)/VLT + k26*tT4LT*fuT4LT - k35*tT3LT*fuT3LT;

%T4LT
dydt(80) = (k25*fT4L - k30*T4LT*fuT4LT)/VLT - k26*T4LT*fuT4LT - k34*T4LT*fuT4LT;

%T3LT
dydt(81) = (k27*fT3L - k31*T3LT*fuT3LT)/VLT + k26*T4LT*fuT4LT - k35*T3LT*fuT3LT;


%Binding Between X and TTR
%XB (Body Blood Compartment)
dydt(82) = k36/VB - k37*XB - k38*XB*TTRB - k42*XB*TBGB + k39*XTTRB + k43*XTBGB + (XT*QT + XRB*QRB + XL*QL - XB*QC)/VB;

%XTTRB (Body Blood Compartment)
dydt(83) = k38*XB*TTRB - k39*XTTRB + (XTTRT*QT + XTTRRB*QRB + XTTRL*QL - XTTRB*QC)/VB;

%XT (Thyroid blood Compartment)
dydt(84) = -k38*XT*TTRT - k42*XT*TBGT + k39*XTTRT + k43*XTBGT + (XB-XT)*QT/VTB;

%XTTRT (Thyroid blood Compartment)
dydt(85) = k38*XT*TTRT - k39*XTTRT + (XTTRB-XTTRT)*QT/VTB;

%XRB (RB blood Compartment)
dydt(86) = -k38*XRB*TTRRB - k42*XRB*TBGRB + k39*XTTRRB + k43*XTBGRB + (XB-XRB)*QRB/VRBB;

%XTTRRB (RB blood Compartment)
dydt(87) = k38*XRB*TTRRB - k39*XTTRRB + (XTTRB-XTTRRB)*QRB/VRBB;

%XL (Liver blood Compartment)
dydt(88) = -k38*XL*TTRL - k42*XL*TBGL + k39*XTTRL + k43*XTBGL + (XB-XL)*QL/VLB;

%XTTRL (Liver blood Compartment)
dydt(89) = k38*XL*TTRL - k39*XTTRL + (XTTRB-XTTRL)*QL/VLB;

%Binding Between X and TBG
%XTBGB (Body Blood Compartment)
dydt(90) = k42*XB*TBGB - k43*XTBGB + (XTBGT*QT + XTBGRB*QRB + XTBGL*QL - XTBGB*QC)/VB;

%XTBGT (Thyroid blood Compartment)
dydt(91) = k42*XT*TBGT - k43*XTBGT + (XTBGB-XTBGT)*QT/VTB;

%XTBGRB (RB blood Compartment)
dydt(92) = k42*XRB*TBGRB - k43*XTBGRB + (XTBGB-XTBGRB)*QRB/VRBB;

%XTBGL (Liver blood Compartment)
dydt(93) = k42*XL*TBGL - k43*XTBGL + (XTBGB-XTBGL)*QL/VLB;

end