clear;
close all ;
clc ;
% warning off;

%% ----------- Initialization ------------------

sign = 1; %neg = 1, pos = 2

fault_duration = 10*60; %Seconds
dt = 1/3600;   % hours
t = [0:dt*60:fault_duration/60];

% Neg
P_ev_initial = -2.5;   %kW
P_bat_initial = -0.5; %kW
P_Load_initial = 1.7;     %kW
SOC_initial_ev = 10;   %kWh    EV Battery SoC at the start of algorithm
SOC_initial_bat = 3;   %kWh     Battery SoC at the start of algorithm
weight_ev = 1;
weight_bat = 1;
weight_load = 0.00001;


% %Pos
% P_ev_initial = -1.8;   %kW
% P_bat_initial = -0.5;  %kW
% P_Load_initial = 2.1;  %kW
% SOC_initial_ev = 31;   %kWh    EV Battery SoC at the start of algorithm
% SOC_initial_bat = 9;   %kWh    Battery SoC at the start of algorithm
% weight_ev = 2.9;
% weight_bat = 1;
% weight_load = 0.00001;

%Flexibility weights


% EV Battery
P_ev_max = 5;          %kW     Max charging rate of EV battery
P_ev_min = -5;         %kW     Max discharging rate of EV battery
SOC_ev_max = 31.5;     %kWh    Max state of charge of EV battery
SOC_ev_min = 3.5;      %kWh    Min state of charge of EV battery


% Battery
P_bat_max = 5;          %kW     Max charging rate of battery
P_bat_min = -5;         %kW     Max discharging rate of battery
SOC_bat_max = 13.5;     %kWh    Max state of charge of battery
SOC_bat_min = 1.5;      %kWh    Min state of charge of battery

% Load load

load_shifting_max_percentage = 0.3; % % of load that can be reduced/increased
P_Load_max = P_Load_initial + P_Load_initial * load_shifting_max_percentage;     %kW
P_Load_min = P_Load_initial - P_Load_initial * load_shifting_max_percentage;     %kW




DP_max = 10;
DP_min = -10;

load('mat_files/f.mat')
% f_test = downsample(f,30);
% f = f_test(1:55);
% t_old_f = linspace(0,fault_duration/60,55);

% if sign == 1
%     
%     f = interp1(t_old_f,f,t,'linear');
%     f = f-50;f = f*0.23;f = f+50; %Negative frequency
% else
%     
%     f = interp1(t_old_f,f,t,'linear')';
%     f = 50-f;f = f*0.25;f = f+50; %Positive frequency
% end

for i = 1:size(t,2)
    DP_dem(i) = Frequency_Controller(f(i),DP_max,DP_min);
end



SOC_ev = SOC_initial_ev;
SOC_bat = SOC_initial_bat;


for i = 1:size(t,2)
    
    A = [dt 0 0;  ...
        -dt 0 0;  ...
        0 dt 0;   ...
        0 -dt 0];
    b = [SOC_ev_max - SOC_ev - P_ev_initial*dt;     ...
        - SOC_ev_min + SOC_ev + P_ev_initial*dt;    ...
        SOC_bat_max - SOC_bat - P_bat_initial*dt;   ...
        - SOC_bat_min + SOC_bat + P_bat_initial*dt];
    
    
    
    if DP_dem(i) >= 0
        
        ub = [P_ev_max - P_ev_initial;  ...
            P_bat_max - P_bat_initial;  ...
            P_Load_max - P_Load_initial];
        
        lb = [0; 0; 0];
        
        
    else
        
        
        ub = [0; 0; 0];
        
        lb = [P_ev_min - P_ev_initial;  ...
            P_bat_min - P_bat_initial; ...
            P_Load_min - P_Load_initial];
        
        
    end
    
    
%     ub = [P_ev_max - P_ev_initial;  ...
%         P_bat_max - P_bat_initial;  ...
%         P_Load_max - P_Load_initial];
%     
%     lb = [P_ev_min - P_ev_initial;  ...
%         P_bat_min - P_bat_initial; ...
%         P_Load_min - P_Load_initial];
    
    
    Aeq = [1 1 1];
    beq = DP_dem(i);
    x0 = [0 0 0];
    
    if DP_dem(i) > 0
        objective = @(x) ( ...
            - weight_ev * ((SOC_ev_max - SOC_ev - (P_ev_initial + x(1) )*dt ) / (SOC_ev_max - SOC_ev_min))        ...
            - weight_bat * ((SOC_bat_max - SOC_bat - (P_bat_initial + x(2) )*dt ) / (SOC_bat_max - SOC_bat_min))   ...
            - weight_load * ((P_Load_max - P_Load_initial - x(3))/(P_Load_max - P_Load_min))                          ...
            );
    elseif DP_dem(i) < 0
        
        objective = @(x) ( ...
            - weight_ev * ((SOC_ev + (P_ev_initial + x(1))*dt - SOC_ev_min) / (SOC_ev_max - SOC_ev_min))        ...
            - weight_bat * ((SOC_bat + (P_bat_initial + x(2))*dt - SOC_bat_min) / (SOC_bat_max - SOC_bat_min))   ...
            - weight_load * ((P_Load_initial + x(3) - P_Load_min)/(P_Load_max - P_Load_min))                          ...
            );
    else
        objective = @(x) (x(1)+ x(2)+ x(3) );
    end
    
    options = optimoptions('fmincon','Display','off','Algorithm','sqp');
    if DP_dem(i) ~= 0
        
        x = fmincon(objective,x0,A,b,Aeq,beq,lb,ub,[],options);
    else
        x(1) = 0;
        x(2) = 0;
        x(3) = 0;
    end
    
    if DP_dem(i) > 0
        flex_Ev(i) = (SOC_ev_max - SOC_ev - (P_ev_initial + x(1) )*dt ) / (SOC_ev_max - SOC_ev_min);
        flex_Bat(i) = (SOC_bat_max - SOC_bat - (P_bat_initial + x(2) )*dt ) / (SOC_bat_max - SOC_bat_min);
        flex_Load(i) = (P_Load_max - P_Load_initial - x(3))/(P_Load_max - P_Load_min);
    elseif DP_dem(i) < 0
        
        flex_Ev(i) =  ((SOC_ev + (P_ev_initial + x(1))*dt - SOC_ev_min) / (SOC_ev_max - SOC_ev_min));
        flex_Bat(i) = ((SOC_bat + (P_bat_initial + x(2))*dt - SOC_bat_min) / (SOC_bat_max - SOC_bat_min));
        flex_Load(i) =  ((P_Load_initial + x(3) - P_Load_min)/(P_Load_max - P_Load_min));
    else
       flex_Ev(i) =  0;
        flex_Bat(i) = 0;
        flex_Load(i) =  0;
    end
    
    DPev(i) = x(1);
    DPbat(i) = x(2);
    DPLoad(i) = x(3);
    Pev(i) = DPev(i) + P_ev_initial;
    Pbat(i) = DPbat(i) + P_bat_initial;
    PLoad(i) = DPLoad(i) + P_Load_initial;
    
    SOC_ev = SOC_ev + Pev(i)*dt;
    SOC_bat = SOC_bat + Pbat(i)*dt;
    
    
    SOC_ev_save(i) = SOC_ev;
    SOC_bat_save(i) = SOC_bat;
    
end


figure('Position',[400 150 650 180])
plot(t,f,'r','LineWidth', 1.5)
% axis([0 fault_duration 49.25 50.25])
ylabel('f_G_R_I_D (Hz)')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.10 .10 .10], 'YColor', [.10 .10 .10], 'YTick', 49.2:0.2:50.7,'LineWidth', 1)
set(gca, 'FontName', 'Helvetica')
xlabel('t (minutes)')


figure('Position',[400 150 650 180])
plot(t,DP_dem,'b','LineWidth', 1.5)
ylabel('ΔP_G_R_I_D (kW)')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.10 .10 .10], 'YColor', [.10 .10 .10], 'YTick', -10:1:10,'LineWidth', 1)
set(gca, 'FontName', 'Helvetica')
xlabel('t (minutes)')


limit_up_ev = (P_ev_max(1))*ones(1,size(t,2));
limit_down_ev = (P_ev_min(1))*ones(1,size(t,2));
limit_up_bat = (P_bat_max(1))*ones(1,size(t,2));
limit_down_bat = (P_bat_min(1))*ones(1,size(t,2));
limit_up_load = (P_Load_max(1))*ones(1,size(t,2));
limit_down_load = (P_Load_min(1))*ones(1,size(t,2));


% figure('Position',[800 150 650 480])
figure
subplot(3,1,1)
plot(t,movmean(Pev,10),'m','LineWidth', 1.5)
title('Power rating')

hold on
ylabel('P_E_V (kW)')
if sign == 1
    plot(t,limit_down_ev,'r--','LineWidth',0.5);
else
    plot(t,limit_up_ev,'r--','LineWidth',0.5);
end
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'on', ...
    'XColor', [.10 .10 .10], 'YColor', [.10 .10 .10], 'YTick', -5:1:5,'LineWidth', 1)
set(gca, 'FontName', 'Helvetica')
legend('P_E_V')

subplot(3,1,2)
plot(t,movmean(Pbat,10),'b','LineWidth', 1.5)
hold on;
ylabel('P_B_A_T (kW)')
if sign == 1
    plot(t,limit_down_bat,'r--','LineWidth',0.5);
else
    plot(t,limit_up_bat,'r--','LineWidth',0.5);
end
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'on', ...
    'XColor', [.10 .10 .10], 'YColor', [.10 .10 .10], 'YTick', -5:1:5,'LineWidth', 1)
set(gca, 'FontName', 'Helvetica')
legend('P_B_A_T')


subplot(3,1,3)
plot(t,PLoad,'g','LineWidth', 1.5)
ylabel('P_L_O_A_D (kW)')
xlabel('t (minutes)')
% axis([0 20 0 3])
hold on
if sign == 1
    plot(t,limit_down_load,'r--','LineWidth',0.5);
else
    plot(t,limit_up_load,'r--','LineWidth',0.5);
end
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'on', ...
    'XColor', [.10 .10 .10], 'YColor', [.10 .10 .10], 'YTick', -5:1:5,'LineWidth', 1)
set(gca, 'FontName', 'Helvetica')
legend('P_L_O_A_D')



limit_up_ev = (SOC_ev_max(1))*ones(1,size(t,2));
limit_down_ev = (SOC_ev_min(1))*ones(1,size(t,2));
limit_up_bat = (SOC_bat_max(1))*ones(1,size(t,2));
limit_down_bat = (SOC_bat_min(1))*ones(1,size(t,2));

figure('Position',[400 150 550 300])

subplot(2,1,1)
plot(t,SOC_ev_save,'m','LineWidth', 1.5)
hold on
if sign == 1
plot(t,limit_down_ev,'r--','LineWidth',0.5);
else
plot(t,limit_up_ev,'r--','LineWidth',0.5);
end
ylabel('SoC_E_V (kWh)')
 legend('SoC_E_V')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'on', ...
    'XColor', [.10 .10 .10], 'YColor', [.10 .10 .10], 'YTick', 0:5:35,'LineWidth', 1)
set(gca, 'FontName', 'Helvetica')
axis([0 10 2 11])
% axis([0 20 29 33])


subplot(2,1,2)
plot(t,SOC_bat_save,'b','LineWidth', 1.5)
hold on
if sign == 1
plot(t,limit_down_bat,'r--','LineWidth',0.5);
else
plot(t,limit_up_bat,'r--','LineWidth',0.5);
end
ylabel('SoC_B_A_T (kWh)')
xlabel('t (minutes)')
legend('SoC_B_A_T')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'on', ...
    'XColor', [.10 .10 .10], 'YColor', [.10 .10 .10], 'YTick', 1:2:16,'LineWidth', 1)
set(gca, 'FontName', 'Helvetica')
axis([0 10 0 4])
% axis([0 20 8.5 14])



% figure
% plot(t, flex_Ev)
% hold on
% plot(t, flex_Bat)
% plot(t, flex_Load)
% legend('flex ev', 'flex bat','flex Load');
% title('Flexibilities')
% axis([0 t(end) 0 1])
% 
% 
% figure
% plot(t, weight_ev * flex_Ev)
% hold on
% plot(t, weight_bat * flex_Bat)
% plot(t, weight_load * flex_Load)
% legend('flex ev', 'flex bat','flex Load');
% title('Weighted Flexibilities')


