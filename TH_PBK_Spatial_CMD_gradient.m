clear all
clc
%--------- Time unit: second, concentration unit: pM, volume unit: L-------%

%% -------------------------- Default Parameters ---------------------------------- %%
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
VPtot = VPC*BW;%Total plasma volume

param.QC   = 15*BW^0.74/3600*(1-HCT);  %Plasma cardiac output
param.QT   = QTC*param.QC;  %Plasma flow rate to Thyroid
param.QRB  = QRBC*param.QC; %Plasma flow rate to RB
param.QL   = QLC*param.QC;  %Plasma flow rate to Liver
param.VTB  = VT*VTBC*(1-HCT);   %Plasma volume in Thyroid
param.VRBB = VRB*VRBBC*(1-HCT); %Plasma volume in RB
param.VRBT = VRB*(1-VRBBC); %Tissue volume in RB
param.VLB  = VL*VLBC*(1-HCT);   %Plasma volume in Liver
param.VLT  = VL*(1-VLBC);   %Tissue volume in Liver
param.VB   = VPtot - param.VTB - param.VRBB - param.VLB;    %Plasma volume of Blood Blood compartment

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
param.k22 = 1.6/14; %pMol/S
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

param.M = 200; %Number of tissue segments

default_param = param;

%% ------------------------- Varrying Parameter Gradient ----------------------------- %% 
% 

arrangement_of_values = 1:1:param.M; 

parameter_indices = 1:14;
parameters = ["k25" "k30" "k27" "k31" "k26" "k34" "k35" "k21" "k28" "k23" "k29" "k24" "k32" "k33"];
parameter_dictionary = containers.Map(parameter_indices, parameters);

for param_index = 1:14
    param_index
    for arrangement = 1:2

        tic

        arrangement
        
        param.k25 = arrayfun(@(current) get_gradient_value(param_index, 1, default_param.k25, current, param.M), arrangement_of_values); % influx_T4_liver_gradients 
        
        param.k30 = arrayfun(@(current) get_gradient_value(param_index, 2, default_param.k30, current, param.M), arrangement_of_values); % efflux_T4_liver_gradients
        
        param.k27 = arrayfun(@(current) get_gradient_value(param_index, 3, default_param.k27, current, param.M), arrangement_of_values); % influx_T3_liver_gradients
        
        param.k31 = arrayfun(@(current) get_gradient_value(param_index, 4, default_param.k31, current, param.M), arrangement_of_values); % efflux_T3_liver_gradient 
        
        param.k26 = arrayfun(@(current) get_gradient_value(param_index, 5, default_param.k26, current, param.M), arrangement_of_values); % conversion_liver_gradients
        
        param.k34 = arrayfun(@(current) get_gradient_value(param_index, 6, default_param.k34, current, param.M), arrangement_of_values); % degradation_T4_liver_gradients
        
        param.k35 = arrayfun(@(current) get_gradient_value(param_index, 7, default_param.k35, current, param.M), arrangement_of_values); % degradation_T3_liver_gradients

        param.k21 = arrayfun(@(current) get_gradient_value(param_index, 8, default_param.k21, current, param.M), arrangement_of_values); % influx_T4_RB_gradients

        param.k28 = arrayfun(@(current) get_gradient_value(param_index, 9, default_param.k28, current, param.M), arrangement_of_values); % efflux_T4_RB_gradients

        param.k23 = arrayfun(@(current) get_gradient_value(param_index, 10, default_param.k23, current, param.M), arrangement_of_values); % influx_T3_RB_gradients

        param.k29 = arrayfun(@(current) get_gradient_value(param_index, 11, default_param.k29, current, param.M), arrangement_of_values); % efflux_T3_RB_gradients

        param.k24 = arrayfun(@(current) get_gradient_value(param_index, 12, default_param.k24, current, param.M), arrangement_of_values); % conversion_RB_gradients

        param.k32 = arrayfun(@(current) get_gradient_value(param_index, 13, default_param.k32, current, param.M), arrangement_of_values); % degradation_T4_RB_gradients
        
        param.k33 = arrayfun(@(current) get_gradient_value(param_index, 14, default_param.k33, current, param.M), arrangement_of_values); % degradation_T3_RB_gradients

        %% ------------------------- Parameter Gradient Plot (Panels A in Fig. 11-13, S5-S16) ----------------------------- %%

        figure(1)
        hold on
        if param_index == 1
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k25, 1:1:param.M))
            end
            plot(1:1:param.M, param.k25)
        elseif param_index == 2
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k30, 1:1:param.M))
            end
            plot(1:1:param.M, param.k30)
        elseif param_index == 3
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k27, 1:1:param.M))
            end
            plot(1:1:param.M, param.k27)
        elseif param_index == 4
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k31, 1:1:param.M))
            end
            plot(1:1:param.M, param.k31)
        elseif param_index == 5
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k26, 1:1:param.M))
            end
            plot(1:1:param.M, param.k26)
        elseif param_index == 6
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k34, 1:1:param.M))
            end
            plot(1:1:param.M, param.k34)
        elseif param_index == 7
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k35, 1:1:param.M))
            end
            plot(1:1:param.M, param.k35)
        elseif param_index == 8
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k21, 1:1:param.M))
            end
            plot(1:1:param.M, param.k21)
        elseif param_index == 9
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k28, 1:1:param.M))
            end
            plot(1:1:param.M, param.k28)
        elseif param_index == 10
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k23, 1:1:param.M))
            end
            plot(1:1:param.M, param.k23)
        elseif param_index == 11
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k29, 1:1:param.M))
            end
            plot(1:1:param.M, param.k29)
        elseif param_index == 12
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k24, 1:1:param.M))
            end
            plot(1:1:param.M, param.k24)
        elseif param_index == 13
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k32, 1:1:param.M))
            end
            plot(1:1:param.M, param.k32)
        else
            if arrangement == 1
                plot(1:1:param.M, arrayfun(@(current) default_param.k33, 1:1:param.M))
            end
            plot(1:1:param.M, param.k33)
        end

        if arrangement == 2
            file_path = sprintf("%s/Fig. XA.fig", parameter_dictionary(param_index));
            saveas(gcf, file_path)
            close all
        end

        
        arrangement_of_values = flip(arrangement_of_values);


        %% ------------------------- Initial Condition ----------------------------- %%
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
        
        init.XB      = 0;
        init.XTTRB   = 0;
        init.XT      = 0;
        init.XTTRT   = 0;
        init.XTBGB   = 0;
        init.XTBGT   = 0;
        
        init.fT4RB(1:param.M,1)   = 15;
        init.T4TBGRB(1:param.M,1) = 0;
        init.T4TTRRB(1:param.M,1) = 0;
        init.T4ALBRB(1:param.M,1) = 0;
        init.fT3RB(1:param.M,1)   = 5;
        init.T3TBGRB(1:param.M,1) = 0;
        init.T3TTRRB(1:param.M,1) = 0;
        init.T3ALBRB(1:param.M,1) = 0;
        init.TBGRB(1:param.M,1)   = 0;
        init.TTRRB(1:param.M,1)   = 0;
        init.ALBRB(1:param.M,1)   = 0;
        init.T4RBT(1:param.M,1)   = 0;
        init.T3RBT(1:param.M,1)   = 0;
        init.XRB(1:param.M,1)     = 0;
        init.XTTRRB(1:param.M,1)  = 0;
        init.XTBGRB(1:param.M,1)  = 0;
        
        init.fT4L(1:param.M,1)    = 15;
        init.T4TBGL(1:param.M,1)  = 0;
        init.T4TTRL(1:param.M,1)  = 0;
        init.T4ALBL(1:param.M,1)  = 0;
        init.fT3L(1:param.M,1)    = 5;
        init.T3TBGL(1:param.M,1)  = 0;
        init.T3TTRL(1:param.M,1)  = 0;
        init.T3ALBL(1:param.M,1)  = 0;
        init.TBGL(1:param.M,1)    = 0;
        init.TTRL(1:param.M,1)    = 0;
        init.ALBL(1:param.M,1)    = 0;
        init.T4LT(1:param.M,1)    = 0;
        init.T3LT(1:param.M,1)    = 0;
        init.XL(1:param.M,1)      = 0;
        init.XTTRL(1:param.M,1)   = 0;
        init.XTBGL(1:param.M,1)   = 0;
        
        init_default = init;
        
        %% ------------------------- Run To Steady-State ----------------------- %%
        y0 = cell2mat(struct2cell(init));
        tspan = [0:10000:1000*24*3600];  %Running for 1000 days
        options = odeset('RelTol',1e-8,'AbsTol',1e-6);
        [t,y] = ode15s('TH_PBK_Spatial_ODE_gradient',tspan, y0, options, param); 
        
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
        
        model.XB      = y(:,23);
        model.XTTRB   = y(:,24);
        model.XT      = y(:,25);
        model.XTTRT   = y(:,26);
        model.XTBGB   = y(:,27);
        model.XTBGT   = y(:,28);
        
        N = 29;  %Starting index of repeating variables
        M = param.M;
        model.fT4RB   = y(:,N:N+M-1);
        model.T4TBGRB = y(:,N+M:N+2*M-1);
        model.T4TTRRB = y(:,N+2*M:N+3*M-1);
        model.T4ALBRB = y(:,N+3*M:N+4*M-1);
        model.fT3RB   = y(:,N+4*M:N+5*M-1);
        model.T3TBGRB = y(:,N+5*M:N+6*M-1);
        model.T3TTRRB = y(:,N+6*M:N+7*M-1);
        model.T3ALBRB = y(:,N+7*M:N+8*M-1);
        model.TBGRB   = y(:,N+8*M:N+9*M-1);
        model.TTRRB   = y(:,N+9*M:N+10*M-1);
        model.ALBRB   = y(:,N+10*M:N+11*M-1);
        model.T4RBT   = y(:,N+11*M:N+12*M-1);
        model.T3RBT   = y(:,N+12*M:N+13*M-1);
        model.XRB     = y(:,N+13*M:N+14*M-1);
        model.XTTRRB  = y(:,N+14*M:N+15*M-1);
        model.XTBGRB  = y(:,N+15*M:N+16*M-1);
        
        model.fT4L    = y(:,N+16*M:N+17*M-1);
        model.T4TBGL  = y(:,N+17*M:N+18*M-1);
        model.T4TTRL  = y(:,N+18*M:N+19*M-1);
        model.T4ALBL  = y(:,N+19*M:N+20*M-1);
        model.fT3L    = y(:,N+20*M:N+21*M-1);
        model.T3TBGL  = y(:,N+21*M:N+22*M-1);
        model.T3TTRL  = y(:,N+22*M:N+23*M-1);
        model.T3ALBL  = y(:,N+23*M:N+24*M-1);
        model.TBGL    = y(:,N+24*M:N+25*M-1);
        model.TTRL    = y(:,N+25*M:N+26*M-1);
        model.ALBL    = y(:,N+26*M:N+27*M-1);
        model.T4LT    = y(:,N+27*M:N+28*M-1);
        model.T3LT    = y(:,N+28*M:N+29*M-1);
        model.XL      = y(:,N+29*M:N+30*M-1);
        model.XTTRL   = y(:,N+30*M:N+31*M-1);
        model.XTBGL   = y(:,N+31*M:N+32*M-1);
        
        
        %-----------Obtaining steady state values-----------%
        % Blood Compartment
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
        
        % RB venous blood (last segment of RB blood) 
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
        
        % RB tissue (last segment)
        T4RBT = model.T4RBT(end);
        T3RBT = model.T3RBT(end);
        
        % Liver venous blood (last segment of Liver blood) 
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
        
        %Liver tissue (last segment)
        T4LT = model.T4LT(end);
        T3LT = model.T3LT(end);
        
        % Calculated metrics for T4 and T3 in Body Blood
        boundT4B = T4TBGB + T4TTRB + T4ALBB;
        totalT4B = boundT4B + freeT4B;
        freeT4Bpercentage = freeT4B/(totalT4B)*100;
        T4TBGBpercentage = T4TBGB/totalT4B*100; 
        T4TTRBpercentage = T4TTRB/totalT4B*100;
        T4ALBBpercentage = T4ALBB/totalT4B*100;
        
        TBGT4Bsatpercentage = T4TBGB/(TBGB+T4TBGB+T3TBGB)*100;
        TTRT4Bsatpercentage = T4TTRB/(TTRB+T4TTRB+T3TTRB)*100;
        ALBT4Bsatpercentage = T4ALBB/(ALBB+T4ALBB+T3ALBB)*100;
        
        boundT3B = T3TBGB + T3TTRB + T3ALBB;
        totalT3B = boundT3B + freeT3B;
        freeT3Bpercentage = freeT3B/(totalT3B)*100;
        T3TBGBpercentage = T3TBGB/totalT3B*100;
        T3TTRBpercentage = T3TTRB/totalT3B*100;
        T3ALBBpercentage = T3ALBB/totalT3B*100;
        
        TBGT3Bsatpercentage = T3TBGB/(TBGB+T4TBGB+T3TBGB)*100;
        TTRT3Bsatpercentage = T3TTRB/(TTRB+T4TTRB+T3TTRB)*100;
        ALBT3Bsatpercentage = T3ALBB/(ALBB+T4ALBB+T3ALBB)*100;
        
        TBGBsatpercentage = TBGT4Bsatpercentage + TBGT3Bsatpercentage;
        TTRBsatpercentage = TTRT4Bsatpercentage + TTRT3Bsatpercentage;
        ALBBsatpercentage = ALBT4Bsatpercentage + ALBT3Bsatpercentage;
        
        % Calculated metrics for T4 and T3 in Thyroid blood
        boundT4T = T4TBGT + T4TTRT + T4ALBT;
        totalT4T = boundT4T + freeT4T;
        freeT4Tpercentage = freeT4T/(totalT4T)*100;
        T4TBGTpercentage = T4TBGT/totalT4T*100;
        T4TTRTpercentage = T4TTRT/totalT4T*100;
        T4ALBTpercentage = T4ALBT/totalT4T*100;
        
        TBGT4Tsatpercentage = T4TBGT/(TBGT+T4TBGT+T3TBGT)*100;
        TTRT4Tsatpercentage = T4TTRT/(TTRT+T4TTRT+T3TTRT)*100;
        ALBT4Tsatpercentage = T4ALBT/(ALBT+T4ALBT+T3ALBT)*100;
        
        boundT3T = T3TBGT + T3TTRT + T3ALBT;
        totalT3T = boundT3T + freeT3T;
        freeT3Tpercentage = freeT3T/(totalT3T)*100;
        T3TBGTpercentage = T3TBGT/totalT3T*100;
        T3TTRTpercentage = T3TTRT/totalT3T*100;
        T3ALBTpercentage = T3ALBT/totalT3T*100;
        
        TBGT3Tsatpercentage = T3TBGT/(TBGT+T4TBGT+T3TBGT)*100;
        TTRT3Tsatpercentage = T3TTRT/(TTRT+T4TTRT+T3TTRT)*100;
        ALBT3Tsatpercentage = T3ALBT/(ALBT+T4ALBT+T3ALBT)*100;
        
        TBGTsatpercentage = TBGT4Tsatpercentage + TBGT3Tsatpercentage;
        TTRTsatpercentage = TTRT4Tsatpercentage + TTRT3Tsatpercentage;
        ALBTsatpercentage = ALBT4Tsatpercentage + ALBT3Tsatpercentage;
        
        % Calculated metrics for T4 and T3 in RB venous blood (last segment of RB blood) 
        boundT4RB = T4TBGRB + T4TTRRB + T4ALBRB;
        totalT4RB = boundT4RB + freeT4RB;
        freeT4RBpercentage = freeT4RB/(totalT4RB)*100;
        T4TBGRBpercentage = T4TBGRB/totalT4RB*100;
        T4TTRRBpercentage = T4TTRRB/totalT4RB*100;
        T4ALBRBpercentage = T4ALBRB/totalT4RB*100;
        
        TBGT4RBsatpercentage = T4TBGRB/(TBGRB+T4TBGRB+T3TBGRB)*100;
        TTRT4RBsatpercentage = T4TTRRB/(TTRRB+T4TTRRB+T3TTRRB)*100;
        ALBT4RBsatpercentage = T4ALBRB/(ALBRB+T4ALBRB+T3ALBRB)*100;
        
        boundT3RB = T3TBGRB + T3TTRRB + T3ALBRB;
        totalT3RB = boundT3RB + freeT3RB;
        freeT3RBpercentage = freeT3RB/(totalT3RB)*100;
        T3TBGRBpercentage = T3TBGRB/totalT3RB*100;
        T3TTRRBpercentage = T3TTRRB/totalT3RB*100;
        T3ALBRBpercentage = T3ALBRB/totalT3RB*100;
        
        TBGT3RBsatpercentage = T3TBGRB/(TBGRB+T4TBGRB+T3TBGRB)*100;
        TTRT3RBsatpercentage = T3TTRRB/(TTRRB+T4TTRRB+T3TTRRB)*100;
        ALBT3RBsatpercentage = T3ALBRB/(ALBRB+T4ALBRB+T3ALBRB)*100;
        
        TBGRBsatpercentage = TBGT4RBsatpercentage + TBGT3RBsatpercentage;
        TTRRBsatpercentage = TTRT4RBsatpercentage + TTRT3RBsatpercentage;
        ALBRBsatpercentage = ALBT4RBsatpercentage + ALBT3RBsatpercentage;
        
        % Calculated metrics for T4 and T3 in Liver venous blood (last segment of Liver blood) 
        boundT4L = T4TBGL + T4TTRL + T4ALBL;
        totalT4L = boundT4L + freeT4L;
        freeT4Lpercentage = freeT4L/(totalT4L)*100;
        T4TBGLpercentage = T4TBGL/totalT4L*100;
        T4TTRLpercentage = T4TTRL/totalT4L*100;
        T4ALBLpercentage = T4ALBL/totalT4L*100;
        
        TBGT4Lsatpercentage = T4TBGL/(TBGL+T4TBGL+T3TBGL)*100;
        TTRT4Lsatpercentage = T4TTRL/(TTRL+T4TTRL+T3TTRL)*100;
        ALBT4Lsatpercentage = T4ALBL/(ALBL+T4ALBL+T3ALBL)*100;
        
        boundT3L = T3TBGL + T3TTRL + T3ALBL;
        totalT3L = boundT3L + freeT3L;
        freeT3Lpercentage = freeT3L/(totalT3L)*100;
        T3TBGLpercentage = T3TBGL/totalT3L*100;
        T3TTRLpercentage = T3TTRL/totalT3L*100;
        T3ALBLpercentage = T3ALBL/totalT3L*100;
        
        TBGT3Lsatpercentage = T3TBGL/(TBGL+T4TBGL+T3TBGL)*100;
        TTRT3Lsatpercentage = T3TTRL/(TTRL+T4TTRL+T3TTRL)*100;
        ALBT3Lsatpercentage = T3ALBL/(ALBL+T4ALBL+T3ALBL)*100;
        
        TBGLsatpercentage = TBGT4Lsatpercentage + TBGT3Lsatpercentage;
        TTRLsatpercentage = TTRT4Lsatpercentage + TTRT3Lsatpercentage;
        ALBLsatpercentage = ALBT4Lsatpercentage + ALBT3Lsatpercentage;
        
        
        % T4 and T3 differnetial between arterial and venous blood concentrations (Table 8 and Fig. 5. Fig.5 was generated in MS Excel using these values)
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
        
        
        % RB blood - all segments
        freeT4RB = model.fT4RB(end, :);
        T4TBGRB = model.T4TBGRB(end, :);
        T4TTRRB = model.T4TTRRB(end, :);
        T4ALBRB = model.T4ALBRB(end, :);
        
        freeT3RB = model.fT3RB(end, :);
        T3TBGRB  = model.T3TBGRB(end, :);
        T3TTRRB  =  model.T3TTRRB(end, :);
        T3ALBRB  = model.T3ALBRB(end, :);
        
        TBGRB = model.TBGRB(end, :);
        TTRRB = model.TTRRB(end, :);
        ALBRB = model.ALBRB(end, :);
        
        % RB tissue - all segments
        T4RBT = model.T4RBT(end, :);
        T3RBT = model.T3RBT(end, :);
        
        % Liver blood - all segments
        freeT4L = model.fT4L(end, :);
        T4TBGL = model.T4TBGL(end, :);
        T4TTRL = model.T4TTRL(end, :);
        T4ALBL = model.T4ALBL(end, :);
        
        freeT3L = model.fT3L(end, :);
        T3TBGL = model.T3TBGL(end, :);
        T3TTRL = model.T3TTRL(end, :);
        T3ALBL = model.T3ALBL(end, :);
        
        TBGL = model.TBGL(end, :);
        TTRL = model.TTRL(end, :);
        ALBL = model.ALBL(end, :);
        
        % Liver tissue - all segments
        T4LT = model.T4LT(end, :);
        T3LT = model.T3LT(end, :);
        
        
        % Body Blood Binding Rates
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
        
        % Thyroid blood Binding Rates
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
        
        % Liver blood Binding Rates
        liver_T4TBG_assoc_rate = param.k1.*freeT4L.*TBGL;
        liver_T4TBG_disassoc_rate = param.k2.*T4TBGL;
        liver_difference_T4TBG_rate = liver_T4TBG_disassoc_rate - liver_T4TBG_assoc_rate; 
        liver_percent_difference_T4TBG_rate = liver_difference_T4TBG_rate./liver_T4TBG_assoc_rate*100;
        
        liver_T4TTR_assoc_rate = param.k3.*freeT4L.*TTRL;
        liver_T4TTR_disassoc_rate = param.k4.*T4TTRL;
        liver_difference_T4TTR_rate = liver_T4TTR_disassoc_rate - liver_T4TTR_assoc_rate; 
        liver_percent_difference_T4TTR_rate = liver_difference_T4TTR_rate./liver_T4TTR_assoc_rate*100;
        
        liver_T4ALB_assoc_rate = param.k5.*freeT4L.*ALBL;
        liver_T4ALB_disassoc_rate = param.k6.*T4ALBL;
        liver_difference_T4ALB_rate = liver_T4ALB_disassoc_rate - liver_T4ALB_assoc_rate; 
        liver_percent_difference_T4ALB_rate = liver_difference_T4ALB_rate./liver_T4ALB_assoc_rate*100;
        
        liver_T3TBG_assoc_rate = param.k7.*freeT3L.*TBGL;
        liver_T3TBG_disassoc_rate = param.k8.*T3TBGL;
        liver_difference_T3TBG_rate = liver_T3TBG_disassoc_rate - liver_T3TBG_assoc_rate; 
        liver_percent_difference_T3TBG_rate = liver_difference_T3TBG_rate./liver_T3TBG_assoc_rate*100;
        
        liver_T3TTR_assoc_rate = param.k9.*freeT3L.*TTRL;
        liver_T3TTR_disassoc_rate = param.k10.*T3TTRL;
        liver_difference_T3TTR_rate = liver_T3TTR_disassoc_rate - liver_T3TTR_assoc_rate; 
        liver_percent_difference_T3TTR_rate = liver_difference_T3TTR_rate./liver_T3TTR_assoc_rate*100;
        
        liver_T3ALB_assoc_rate = param.k11.*freeT3L.*ALBL;
        liver_T3ALB_disassoc_rate = param.k12.*T3ALBL;
        liver_difference_T3ALB_rate = liver_T3ALB_disassoc_rate - liver_T3ALB_assoc_rate; 
        liver_percent_difference_T3ALB_rate = liver_difference_T3ALB_rate./liver_T3ALB_assoc_rate*100;
        
        % RB blood Binding Rates
        RB_T4TBG_assoc_rate = param.k1.*freeT4RB.*TBGRB;
        RB_T4TBG_disassoc_rate = param.k2.*T4TBGRB;
        RB_difference_T4TBG_rate = RB_T4TBG_disassoc_rate - RB_T4TBG_assoc_rate; 
        RB_percent_difference_T4TBG_rate = RB_difference_T4TBG_rate./RB_T4TBG_assoc_rate*100;
        
        RB_T4TTR_assoc_rate = param.k3.*freeT4RB.*TTRRB;
        RB_T4TTR_disassoc_rate = param.k4.*T4TTRRB;
        RB_difference_T4TTR_rate = RB_T4TTR_disassoc_rate - RB_T4TTR_assoc_rate; 
        RB_percent_difference_T4TTR_rate = RB_difference_T4TTR_rate./RB_T4TTR_assoc_rate*100;
        
        RB_T4ALB_assoc_rate = param.k5.*freeT4RB.*ALBRB;
        RB_T4ALB_disassoc_rate = param.k6.*T4ALBRB;
        RB_difference_T4ALB_rate = RB_T4ALB_disassoc_rate - RB_T4ALB_assoc_rate; 
        RB_percent_difference_T4ALB_rate = RB_difference_T4ALB_rate./RB_T4ALB_assoc_rate*100;
        
        RB_T3TBG_assoc_rate = param.k7.*freeT3RB.*TBGRB;
        RB_T3TBG_disassoc_rate = param.k8.*T3TBGRB;
        RB_difference_T3TBG_rate = RB_T3TBG_disassoc_rate - RB_T3TBG_assoc_rate; 
        RB_percent_difference_T3TBG_rate = RB_difference_T3TBG_rate./RB_T3TBG_assoc_rate*100;
        
        RB_T3TTR_assoc_rate = param.k9.*freeT3RB.*TTRRB;
        RB_T3TTR_disassoc_rate = param.k10.*T3TTRRB;
        RB_difference_T3TTR_rate = RB_T3TTR_disassoc_rate - RB_T3TTR_assoc_rate; 
        RB_percent_difference_T3TTR_rate = RB_difference_T3TTR_rate./RB_T3TTR_assoc_rate*100;
        
        RB_T3ALB_assoc_rate = param.k11.*freeT3RB.*ALBRB;
        RB_T3ALB_disassoc_rate = param.k12.*T3ALBRB;
        RB_difference_T3ALB_rate = RB_T3ALB_disassoc_rate - RB_T3ALB_assoc_rate; 
        RB_percent_difference_T3ALB_rate = RB_difference_T3ALB_rate./RB_T3ALB_assoc_rate*100;
        
        
        %% ------------------------- Plotting Results ---------------------------------- %%

        %Free T4 in Liver and RB blood (Panels B in Fig. 11, 13, S5, S7-S16, and panel D in Fig. S6)
        figure(30)
        hold on
        plot([0:1:param.M],[model.fT4B(end),model.fT4L(end,:)])
        plot([0:1:param.M],[model.fT4B(end),model.fT4RB(end,:)])
        file_path = sprintf("%s %d/Fig. XB.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %T4TBG in Liver and RB blood (not used)
        figure(31)
        hold on
        plot([0:1:param.M],[model.T4TBGB(end),model.T4TBGL(end,:)])
        plot([0:1:param.M],[model.T4TBGB(end),model.T4TBGRB(end,:)])
        file_path = sprintf("%s %d/Fig. X.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %T4TTR in Liver and RB blood (not used)
        figure(32)
        hold on
        plot([0:1:param.M],[model.T4TTRB(end),model.T4TTRL(end,:)])
        plot([0:1:param.M],[model.T4TTRB(end),model.T4TTRRB(end,:)])
        file_path = sprintf("%s %d/Fig. X.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %T4ALB in Liver and RB blood (not used)
        figure(33)
        hold on
        plot([0:1:param.M],[model.T4ALBB(end),model.T4ALBL(end,:)])
        plot([0:1:param.M],[model.T4ALBB(end),model.T4ALBRB(end,:)])
        file_path = sprintf("%s %d/Fig. X.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %Total T4 in Liver and RB blood (Panels C in Fig. 11, 13, S5, S7-S16, and panel E in Fig. S6)
        figure(34)
        hold on
        plot([0:1:param.M],[model.fT4B(end)+model.T4TBGB(end)+model.T4TTRB(end)+model.T4ALBB(end), model.fT4L(end,:)+model.T4TBGL(end,:)+model.T4TTRL(end,:)+model.T4ALBL(end,:)])
        plot([0:1:param.M],[model.fT4B(end)+model.T4TBGB(end)+model.T4TTRB(end)+model.T4ALBB(end), model.fT4RB(end,:)+model.T4TBGRB(end,:)+model.T4TTRRB(end,:)+model.T4ALBRB(end,:)])
        file_path = sprintf("%s %d/Fig. XC.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        
        %Free T3 in Liver and RB blood (Panels G in Fig. 11, 13, S5, S7-S16, and panel D in Fig. 12)
        figure(35)
        hold on
        plot([0:1:param.M],[model.fT3B(end),model.fT3L(end,:)])
        plot([0:1:param.M],[model.fT3B(end),model.fT3RB(end,:)])
        file_path = sprintf("%s %d/Fig. XG.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %T3TBG in Liver and RB blood (not used)
        figure(36)
        hold on
        plot([0:1:param.M],[model.T3TBGB(end),model.T3TBGL(end,:)])
        plot([0:1:param.M],[model.T3TBGB(end),model.T3TBGRB(end,:)])
        file_path = sprintf("%s %d/Fig. X.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %T3TTR in Liver and RB blood (not used)
        figure(37)
        hold on
        plot([0:1:param.M],[model.T3TTRB(end),model.T3TTRL(end,:)])
        plot([0:1:param.M],[model.T3TTRB(end),model.T3TTRRB(end,:)])
        file_path = sprintf("%s %d/Fig. X.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %T3ALB in Liver and RB blood (not used)
        figure(38)
        hold on
        plot([0:1:param.M],[model.T3ALBB(end),model.T3ALBL(end,:)])
        plot([0:1:param.M],[model.T3ALBB(end),model.T3ALBRB(end,:)])
        file_path = sprintf("%s %d/Fig. X.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %Total T3 in Liver and RB blood (Panels H in Fig. 11, 13, S5, S7-S16, and panel E in Fig. 12)
        figure(39)
        hold on
        plot([0:1:param.M],[model.fT3B(end)+model.T3TBGB(end)+model.T3TTRB(end)+model.T3ALBB(end), model.fT3L(end,:)+model.T3TBGL(end,:)+model.T3TTRL(end,:)+model.T3ALBL(end,:)])
        plot([0:1:param.M],[model.fT3B(end)+model.T3TBGB(end)+model.T3TTRB(end)+model.T3ALBB(end), model.fT3RB(end,:)+model.T3TBGRB(end,:)+model.T3TTRRB(end,:)+model.T3ALBRB(end,:)])
        file_path = sprintf("%s %d/Fig. XH.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %T4 in Liver and RB tissue (Panels D in Fig. 11, 13, S5, S7-S16, and panel F in Fig. S6)
        figure(301)
        hold on
        plot([1:1:param.M],model.T4LT(end,:))
        plot([1:1:param.M],model.T4RBT(end,:))
        file_path = sprintf("%s %d/Fig. XD.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %T3 in Liver and RB tissue  (Panels I in Fig. 11, 13, S5, S7-S16, and panel F in Fig. 12)
        figure(302)
        hold on
        plot([1:1:param.M],model.T3LT(end,:))
        plot([1:1:param.M],model.T3RBT(end,:))
        file_path = sprintf("%s %d/Fig. XI.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %Free TBG in Liver and RB Blood (not used)
        figure(303)
        hold on
        plot([0:1:param.M],[model.TBGB(end),model.TBGL(end,:)])
        plot([0:1:param.M],[model.TBGB(end),model.TBGRB(end,:)])
        file_path = sprintf("%s %d/Fig. X.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %Free TTR in Liver and RB Blood (not used)
        figure(304)
        hold on
        plot([0:1:param.M],[model.TTRB(end),model.TTRL(end,:)])
        plot([0:1:param.M],[model.TTRB(end),model.TTRRB(end,:)])
        file_path = sprintf("%s %d/Fig. X.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %Free ALB in Liver and RB Blood (not used)
        figure(305)
        hold on
        plot([0:1:param.M],[model.ALBB(end),model.ALBL(end,:)])
        plot([0:1:param.M],[model.ALBB(end),model.ALBRB(end,:)])
        file_path = sprintf("%s %d/Fig. X.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
               
        %Absolute Drop of T4 in Liver blood (Panels E and F in Fig. 11, 13, S5, S7-S9, and panels B and C in Fig. S6)
        figure(41)
        hold on
        plot([1:1:param.M], [model.T4ALBB(end)-model.T4ALBL(end,1) + model.T4TTRB(end)-model.T4TTRL(end,1) + model.T4TBGB(end)-model.T4TBGL(end,1) + model.fT4B(end)-model.fT4L(end,1), - model.T4TBGL(end,2:param.M) + model.T4TBGL(end,1:param.M-1) - model.T4ALBL(end,2:param.M) + model.T4ALBL(end,1:param.M-1) - model.T4TTRL(end,2:param.M) + model.T4TTRL(end,1:param.M-1) - model.fT4L(end,2:param.M) + model.fT4L(end,1:param.M-1)])
        plot([1:1:param.M], [model.T4ALBB(end)-model.T4ALBL(end,1), - model.T4ALBL(end,2:param.M) + model.T4ALBL(end,1:param.M-1)])
        plot([1:1:param.M], [model.T4TTRB(end)-model.T4TTRL(end,1), - model.T4TTRL(end,2:param.M) + model.T4TTRL(end,1:param.M-1)])
        plot([1:1:param.M], [model.T4TBGB(end)-model.T4TBGL(end,1), - model.T4TBGL(end,2:param.M) + model.T4TBGL(end,1:param.M-1)])
        plot([1:1:param.M], [model.fT4B(end)-model.fT4L(end,1), - model.fT4L(end,2:param.M) + model.fT4L(end,1:param.M-1)])
        file_path = sprintf("%s %d/Fig. XEF.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %Absolute Drop of T4 in RB blood (Panels E and F in Fig. S10-S16) 
        figure(42)
        hold on
        plot([1:1:param.M], [model.T4ALBB(end)-model.T4ALBRB(end,1) + model.T4TTRB(end)-model.T4TTRRB(end,1) + model.T4TBGB(end)-model.T4TBGRB(end,1) + model.fT4B(end)-model.fT4RB(end,1), - model.T4TBGRB(end,2:param.M) + model.T4TBGRB(end,1:param.M-1) - model.T4ALBRB(end,2:param.M) + model.T4ALBRB(end,1:param.M-1) - model.T4TTRRB(end,2:param.M) + model.T4TTRRB(end,1:param.M-1) - model.fT4RB(end,2:param.M) + model.fT4RB(end,1:param.M-1)])
        plot([1:1:param.M], [model.T4ALBB(end)-model.T4ALBRB(end,1), - model.T4ALBRB(end,2:param.M) + model.T4ALBRB(end,1:param.M-1)])
        plot([1:1:param.M], [model.T4TTRB(end)-model.T4TTRRB(end,1), - model.T4TTRRB(end,2:param.M) + model.T4TTRRB(end,1:param.M-1)])
        plot([1:1:param.M], [model.T4TBGB(end)-model.T4TBGRB(end,1), - model.T4TBGRB(end,2:param.M) + model.T4TBGRB(end,1:param.M-1)])
        plot([1:1:param.M], [model.fT4B(end)-model.fT4RB(end,1), - model.fT4RB(end,2:param.M) + model.fT4RB(end,1:param.M-1)])
        file_path = sprintf("%s %d/Fig. XEF.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %Absolute Drop of T3 in Liver blood (Panels J and K in Fig. 11, 13, S5, S7-S9, and panels B and C in Fig. 12)
        figure(43)
        hold on
        plot([1:1:param.M], [model.T3ALBB(end)-model.T3ALBL(end,1) + model.T3TTRB(end)-model.T3TTRL(end,1) + model.T3TBGB(end)-model.T3TBGL(end,1) + model.fT3B(end)-model.fT3L(end,1), - model.T3TBGL(end,2:param.M) + model.T3TBGL(end,1:param.M-1) - model.T3ALBL(end,2:param.M) + model.T3ALBL(end,1:param.M-1) - model.T3TTRL(end,2:param.M) + model.T3TTRL(end,1:param.M-1) - model.fT3L(end,2:param.M) + model.fT3L(end,1:param.M-1)])
        plot([1:1:param.M], [model.T3ALBB(end)-model.T3ALBL(end,1), - model.T3ALBL(end,2:param.M) + model.T3ALBL(end,1:param.M-1)])
        plot([1:1:param.M], [model.T3TTRB(end)-model.T3TTRL(end,1), - model.T3TTRL(end,2:param.M) + model.T3TTRL(end,1:param.M-1)])
        plot([1:1:param.M], [model.T3TBGB(end)-model.T3TBGL(end,1), - model.T3TBGL(end,2:param.M) + model.T3TBGL(end,1:param.M-1)])
        plot([1:1:param.M], [model.fT3B(end)-model.fT3L(end,1), - model.fT3L(end,2:param.M) + model.fT3L(end,1:param.M-1)])
        file_path = sprintf("%s %d/Fig. XJK.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        
        % %Absolute Drop of T3 in RB blood (Panels J and K in Fig. S10-S16) 
        figure(44)
        hold on
        plot([1:1:param.M], [model.T3ALBB(end)-model.T3ALBRB(end,1) + model.T3TTRB(end)-model.T3TTRRB(end,1) + model.T3TBGB(end)-model.T3TBGRB(end,1) + model.fT3B(end)-model.fT3RB(end,1), - model.T3TBGRB(end,2:param.M) + model.T3TBGRB(end,1:param.M-1) - model.T3ALBRB(end,2:param.M) + model.T3ALBRB(end,1:param.M-1) - model.T3TTRRB(end,2:param.M) + model.T3TTRRB(end,1:param.M-1) - model.fT3RB(end,2:param.M) + model.fT3RB(end,1:param.M-1)])
        plot([1:1:param.M], [model.T3ALBB(end)-model.T3ALBRB(end,1), - model.T3ALBRB(end,2:param.M) + model.T3ALBRB(end,1:param.M-1)])
        plot([1:1:param.M], [model.T3TTRB(end)-model.T3TTRRB(end,1), - model.T3TTRRB(end,2:param.M) + model.T3TTRRB(end,1:param.M-1)])
        plot([1:1:param.M], [model.T3TBGB(end)-model.T3TBGRB(end,1), - model.T3TBGRB(end,2:param.M) + model.T3TBGRB(end,1:param.M-1)])
        plot([1:1:param.M], [model.fT3B(end)-model.fT3RB(end,1), - model.fT3RB(end,2:param.M) + model.fT3RB(end,1:param.M-1)])
        file_path = sprintf("%s %d/Fig. XJK.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        
        %% ------------------------- Plotting Results (not used) ---------------------------------- %%
        
        %Liver blood T4TBG assocation and dissociation
        figure(51)
        hold on
        yyaxis left
        plot([0:1:param.M],[body_blood_T4TBG_assoc_rate, liver_T4TBG_assoc_rate])
        plot([0:1:param.M], [body_blood_T4TBG_disassoc_rate, liver_T4TBG_disassoc_rate])
        xlabel('Tissue Segment')
        ylabel('Rate (pM/S)')
        yyaxis right
        plot([0:1:param.M], [body_blood_percent_difference_T4TBG_rate, liver_percent_difference_T4TBG_rate])
        ylabel('Percent Difference')
        file_path = sprintf("/Users/maxbagga/Library/CloudStorage/OneDrive-EmoryUniversity/TH binding modeling/TH PBK Models/1_2_PBK/2_PBK Spatial/Intermediate result figures/Gradient Analysis/%s %d/Fig. 5A.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %Liver blood T4TTR assocation and dissociation
        figure(52)
        hold on
        yyaxis left
        plot([0:1:param.M], [body_blood_T4TTR_assoc_rate, liver_T4TTR_assoc_rate])
        plot([0:1:param.M], [body_blood_T4TTR_disassoc_rate, liver_T4TTR_disassoc_rate])
        xlabel('Tissue Segment')
        ylabel('Rate (pM/S)')
        yyaxis right
        plot([0:1:param.M], [body_blood_percent_difference_T4TTR_rate, liver_percent_difference_T4TTR_rate])
        ylabel('Percent Difference')
        file_path = sprintf("/Users/maxbagga/Library/CloudStorage/OneDrive-EmoryUniversity/TH binding modeling/TH PBK Models/1_2_PBK/2_PBK Spatial/Intermediate result figures/Gradient Analysis/%s %d/Fig. 5B.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %Liver blood T4ALB assocation and dissociation
        figure(53)
        hold on
        yyaxis left
        plot([0:1:param.M], [body_blood_T4ALB_assoc_rate, liver_T4ALB_assoc_rate])
        plot([0:1:param.M], [body_blood_T4ALB_disassoc_rate, liver_T4ALB_disassoc_rate])
        xlabel('Tissue Segment')
        ylabel('Rate (pM/S)')
        yyaxis right
        plot([0:1:param.M], [body_blood_percent_difference_T4ALB_rate liver_percent_difference_T4ALB_rate])
        ylabel('Percent Difference')
        file_path = sprintf("/Users/maxbagga/Library/CloudStorage/OneDrive-EmoryUniversity/TH binding modeling/TH PBK Models/1_2_PBK/2_PBK Spatial/Intermediate result figures/Gradient Analysis/%s %d/Fig. 5C.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %RB blood T4TBG assocation and dissociation
        figure(54)
        hold on
        yyaxis left
        plot([0:1:param.M], [body_blood_T4TBG_assoc_rate, RB_T4TBG_assoc_rate])
        plot([0:1:param.M], [body_blood_T4TBG_disassoc_rate, RB_T4TBG_disassoc_rate])
        xlabel('Tissue Segment')
        ylabel('Rate (pM/S)')
        yyaxis right
        plot([0:1:param.M], [body_blood_percent_difference_T4TBG_rate, RB_percent_difference_T4TBG_rate])
        ylabel('Percent Difference')
        file_path = sprintf("/Users/maxbagga/Library/CloudStorage/OneDrive-EmoryUniversity/TH binding modeling/TH PBK Models/1_2_PBK/2_PBK Spatial/Intermediate result figures/Gradient Analysis/%s %d/Fig. 5D.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %RB blood T4TTR assocation and dissociation
        figure(55)
        hold on
        yyaxis left
        plot([0:1:param.M], [body_blood_T4TTR_assoc_rate, RB_T4TTR_assoc_rate])
        plot([0:1:param.M], [body_blood_T4TTR_disassoc_rate, RB_T4TTR_disassoc_rate])
        xlabel('Tissue Segment')
        ylabel('Rate (pM/S)')
        yyaxis right
        plot([0:1:param.M], [body_blood_percent_difference_T4TTR_rate, RB_percent_difference_T4TTR_rate])
        ylabel('Percent Difference')
        file_path = sprintf("/Users/maxbagga/Library/CloudStorage/OneDrive-EmoryUniversity/TH binding modeling/TH PBK Models/1_2_PBK/2_PBK Spatial/Intermediate result figures/Gradient Analysis/%s %d/Fig. 5E.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %RB blood T4ALB assocation and dissociation
        figure(56)
        hold on
        yyaxis left
        plot([0:1:param.M], [body_blood_T4ALB_assoc_rate, RB_T4ALB_assoc_rate])
        plot([0:1:param.M], [body_blood_T4ALB_disassoc_rate, RB_T4ALB_disassoc_rate])
        xlabel('Tissue Segment')
        ylabel('Rate (pM/S)')
        yyaxis right
        plot([0:1:param.M], [body_blood_percent_difference_T4ALB_rate, RB_percent_difference_T4ALB_rate])
        ylabel('Percent Difference')
        file_path = sprintf("/Users/maxbagga/Library/CloudStorage/OneDrive-EmoryUniversity/TH binding modeling/TH PBK Models/1_2_PBK/2_PBK Spatial/Intermediate result figures/Gradient Analysis/%s %d/Fig. 5F.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        
        %Liver blood T3TBG assocation and dissociation
        figure(61)
        hold on
        yyaxis left
        plot([0:1:param.M], [body_blood_T3TBG_assoc_rate, liver_T3TBG_assoc_rate])
        plot([0:1:param.M], [body_blood_T3TBG_disassoc_rate, liver_T3TBG_disassoc_rate])
        xlabel('Tissue Segment')
        ylabel('Rate (pM/S)')
        yyaxis right
        plot([0:1:param.M], [body_blood_percent_difference_T3TBG_rate, liver_percent_difference_T3TBG_rate])
        ylabel('Percent Difference')
        file_path = sprintf("/Users/maxbagga/Library/CloudStorage/OneDrive-EmoryUniversity/TH binding modeling/TH PBK Models/1_2_PBK/2_PBK Spatial/Intermediate result figures/Gradient Analysis/%s %d/Fig. 6A.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %Liver blood T3TTR assocation and dissociation
        figure(62)
        hold on
        yyaxis left
        plot([0:1:param.M], [body_blood_T3TTR_assoc_rate, liver_T3TTR_assoc_rate])
        plot([0:1:param.M], [body_blood_T3TTR_disassoc_rate, liver_T3TTR_disassoc_rate])
        xlabel('Tissue Segment')
        ylabel('Rate (pM/S)')
        yyaxis right
        plot([0:1:param.M], [body_blood_percent_difference_T3TTR_rate, liver_percent_difference_T3TTR_rate])
        ylabel('Percent Difference')
        file_path = sprintf("/Users/maxbagga/Library/CloudStorage/OneDrive-EmoryUniversity/TH binding modeling/TH PBK Models/1_2_PBK/2_PBK Spatial/Intermediate result figures/Gradient Analysis/%s %d/Fig. 6B.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %Liver blood T3ALB assocation and dissociation
        figure(63)
        hold on
        yyaxis left
        plot([0:1:param.M], [body_blood_T3ALB_assoc_rate, liver_T3ALB_assoc_rate])
        plot([0:1:param.M], [body_blood_T3ALB_disassoc_rate, liver_T3ALB_disassoc_rate])
        xlabel('Tissue Segment')
        ylabel('Rate (pM/S)')
        yyaxis right
        plot([0:1:param.M], [body_blood_percent_difference_T3ALB_rate, liver_percent_difference_T3ALB_rate])
        ylabel('Percent Difference')
        file_path = sprintf("/Users/maxbagga/Library/CloudStorage/OneDrive-EmoryUniversity/TH binding modeling/TH PBK Models/1_2_PBK/2_PBK Spatial/Intermediate result figures/Gradient Analysis/%s %d/Fig. 6C.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %RB blood T3TBG assocation and dissociation
        figure(64)
        hold on
        yyaxis left
        plot([0:1:param.M], [body_blood_T3TBG_assoc_rate, RB_T3TBG_assoc_rate])
        plot([0:1:param.M], [body_blood_T3TBG_disassoc_rate, RB_T3TBG_disassoc_rate])
        xlabel('Tissue Segment')
        ylabel('Rate (pM/S)')
        yyaxis right
        plot([0:1:param.M], [body_blood_percent_difference_T3TBG_rate, RB_percent_difference_T3TBG_rate])
        ylabel('Percent Difference')
        file_path = sprintf("/Users/maxbagga/Library/CloudStorage/OneDrive-EmoryUniversity/TH binding modeling/TH PBK Models/1_2_PBK/2_PBK Spatial/Intermediate result figures/Gradient Analysis/%s %d/Fig. 6D.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %RB blood T3TTR assocation and dissociation
        figure(65)
        hold on
        yyaxis left
        plot([0:1:param.M], [body_blood_T3TTR_assoc_rate, RB_T3TTR_assoc_rate])
        plot([0:1:param.M], [body_blood_T3TTR_disassoc_rate, RB_T3TTR_disassoc_rate])
        xlabel('Tissue Segment')
        ylabel('Rate (pM/S)')
        yyaxis right
        plot([0:1:param.M], [body_blood_percent_difference_T3TTR_rate, RB_percent_difference_T3TTR_rate])
        ylabel('Percent Difference')
        file_path = sprintf("/Users/maxbagga/Library/CloudStorage/OneDrive-EmoryUniversity/TH binding modeling/TH PBK Models/1_2_PBK/2_PBK Spatial/Intermediate result figures/Gradient Analysis/%s %d/Fig. 6E.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)
        
        %RB blood T3ALB assocation and dissociation
        figure(66)
        hold on
        yyaxis left
        plot([0:1:param.M], [body_blood_T3ALB_assoc_rate, RB_T3ALB_assoc_rate])
        plot([0:1:param.M], [body_blood_T3ALB_disassoc_rate, RB_T3ALB_disassoc_rate])
        xlabel('Tissue Segment')
        ylabel('Rate (pM/S)')
        yyaxis right
        plot([0:1:param.M], [body_blood_percent_difference_T3ALB_rate, RB_percent_difference_T3ALB_rate])
        ylabel('Percent Difference')
        file_path = sprintf("/Users/maxbagga/Library/CloudStorage/OneDrive-EmoryUniversity/TH binding modeling/TH PBK Models/1_2_PBK/2_PBK Spatial/Intermediate result figures/Gradient Analysis/%s %d/Fig. 6F.fig", parameter_dictionary(param_index), arrangement);
        saveas(gcf, file_path)

        close all

        toc

    end
end

%% ------------------------- Gradient Function ----------------------------- %%

function gradient = get_gradient_value(current_index, index_of_param, ...
    parameter, current, M)
% Calculate gradient for parameter
    % To generate exponential gradient
    % min_factor = 0.1;
    % geo_ratio = 1.01830062850419; % Found by solving exponential sum formula for the geometric ratio where the number of terms is M, sum is target parameter value * M, and initial value is min_factor * parameter
    % current_value = @(parameter, current) parameter*min_factor * (geo_ratio^(current - 1));
    
    % To generate linear gradient
    fold = 10;
    variation_difference = (fold - 1)/(fold+1);

    if index_of_param == current_index
        gradient = parameter * (1 - variation_difference + ...
            (current - 1) * 2 * variation_difference / (M - 1));
    else
        gradient = parameter;
    end
end
