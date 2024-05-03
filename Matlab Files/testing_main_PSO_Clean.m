clear all;
close all;
clc;


dt=0.5;
T=24;
t = 0:dt:T;
T_end = 24;

% CostOfPower(1:14)= 0.07897; CostOfPower(15:46)= 0.11058; CostOfPower(47:49)= 0.07897; %€/kwh
COP_Ther_Var_struct = load('cop_ther_var.mat');
COP_Ther_Var=COP_Ther_Var_struct.vq1;               %Cost of Power - Summer - Variable
CostOfPower = COP_Ther_Var;


numParticles = 30;
maxEpochs = 250;

Pmax_ev = 5;
Pmin_ev = -5;
EV_bat_n = 0.85;
SoC_initial_ev = 19.92;
dt = 0.5;
Ev_bat_capacity = 35;
SoC_target_ev = 20;
SoCmax_ev = 35*0.9;
SoCmin_ev = 35*0.1;

Tmax = 25;
Tmin = 23;
T_in_initial = 25;
V = 320;
C = 1;
rho = 1.2;
Uwall = 0.00048;
Uwindow = 0.0029;
Fwall = 224;
Fwindow = 20;
aw = 18.63;
Rsej = 0.04;
taf_window = 1.1;
SC = 0.54;
QEC_max = 5;
QEC_min = 0;
COP = 3.5;
solar_radiation_struct = load('solar_radiation.mat');
outdoor_temperature_struct = load('summer_daily_temperature.mat');
outdoor_temperature = outdoor_temperature_struct.data(1:49);
solar_radiation = solar_radiation_struct.y(1:49);
heatcool_sign = -1;


Zeta = (rho*C*V)/(dt*3600);     %dt is converted to seconds in order to match the rest of the units
% Alpha = 1/Zeta;
% Beta = (Zeta-(Uwall*Fwall)-(Uwindow*Fwindow))/Zeta;
Qin = 0.2*ones(1,length(t));
Qsw = aw * Rsej * Uwall * Fwall * solar_radiation;
Qsg = 0.001*taf_window * SC * Fwindow * solar_radiation;
% Qin = zeros(1,length(t));
% Qsw = zeros(1,length(t));
% Qsg = zeros(1,length(t));


U = ((Uwall*Fwall+Uwindow*Fwindow)/Zeta)*outdoor_temperature+(Qin+Qsw+Qsg)/Zeta;
U_initial=((Uwall*Fwall+Uwindow*Fwindow)/Zeta)*outdoor_temperature(1);



for i=1: (numParticles*(maxEpochs+1))+1
    cost = psoClean(numParticles,maxEpochs);
end


%% ------------- EV Calculations -------------
PPEV = particlesTable(i_min_cost,76:124);
PPEV_actual = zeros(1,length(PPEV))';
for z=1:length(PPEV_actual)
    if PPEV(z)>=0
        PPEV_actual(z)=PPEV(z)*EV_bat_n;
    else
        PPEV_actual(z)=PPEV(z)/EV_bat_n;
    end
end

state_of_charge_ev = zeros(1,length(PPEV));
state_of_charge_ev(1)=SoC_initial_ev;
state_of_charge_ev(2:end)=SoC_initial_ev+cumsum(PPEV_actual(1:end-1))*dt;

%% ------------- Heat System Calculations -------------


QEC_not_interp = particlesTable(i_min_cost,125:149);  % TODO Fix number of variables once you add more power arrays to calculate
QEC = interp1([0:1:24],QEC_not_interp,[0:dt:24],'spline');
QEC_thermal = QEC*COP;


indoor_temperature=zeros(1,length(QEC_thermal));
indoor_temperature(1) =T_in_initial;

for i=2:length(QEC)
    indoor_temperature(i) = indoor_temperature(i-1) + heatcool_sign*QEC_thermal(i-1)/Zeta + ...
        (Uwall*Fwall*(outdoor_temperature(i-1)-indoor_temperature(i-1))+ ...
        Uwindow*Fwindow*(outdoor_temperature(i-1)-indoor_temperature(i-1))+ ...
        Qin(i-1)+Qsw(i-1)+Qsg(i-1))/Zeta;
end


figure
xBox = [0, 24, 24, 0];
yBox = [0, 0, 5, 5];
patch(xBox, yBox, 'black', 'FaceColor', 'green', 'FaceAlpha', 0.1);
xBox = [0 24 24 0];
yBox = [0 0 -5 -5];
patch(xBox, yBox, 'black', 'FaceColor', 'red', 'FaceAlpha', 0.1);
hold on;
plot(t,PPEV,'b','LineWidth', 1.5);
%                 plot(t,PPEV_actual,'r','LineWidth', 1.5);
axis([0 24 -5 5])
title('Power Exchange of EV battery')
xlabel('t (hours)')
ylabel('P_P_E_V (kW)')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', -5:2.5:5,'LineWidth', 1)
set(gca, 'FontName', 'Helvetica')
xticks(0:4:24)
xticklabels({'00:00','04:00','08:00','12:00','16:00','20:00','24:00'})


figure
plot(t,state_of_charge_ev,'m','LineWidth', 1.5)
hold on;
limitline_top = SoCmax_ev*ones(1,length(state_of_charge_ev));
plot(t,limitline_top,'r--');
limitline_bottom = SoCmin_ev*ones(1,length(state_of_charge_ev));
plot(t,limitline_bottom,'r--');
hold off;
axis([0 24 0 35])
xlabel('t (hours)')
ylabel('E (kWh)')
title('Energy of EV Battery')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', 0:10:Ev_bat_capacity,'LineWidth', 1)
xticks(0:4:24)
xticklabels({'00:00','04:00','08:00','12:00','16:00','20:00','24:00'})



figure
plot(t,CostOfPower,'r','LineWidth', 1.5)
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YTick', min(CostOfPower):0.01:max(CostOfPower),'LineWidth', 1)
set(gca, 'FontName', 'Helvetica')
hold on;
axis([0 24 0.07 0.14])
ylabel('Cost (€ / kWh)')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', 0.07:0.01:0.14,'LineWidth', 1)
title('Cost of Power')
xlabel('t (hours)')
xticks(0:4:24)
xticklabels({'00:00','04:00','08:00','12:00','16:00','20:00','24:00'})



figure
plot(t,QEC,'c','LineWidth', 1.5)
xlabel('t (hours)')
ylabel('P_Q_E_C (kW)')
title('Electric Chiller Power')
axis([0 24 QEC_min QEC_max])
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', QEC_min:1:QEC_max,'LineWidth', 1)
xticks(0:4:24)
xticklabels({'00:00','04:00','08:00','12:00','16:00','20:00','24:00'})

figure
plot(t,outdoor_temperature,'r','LineWidth', 1.5)
hold on;
plot(t,indoor_temperature,'b','LineWidth', 1.5)
limitline_top = Tmax*ones(1,length(t));
plot(t,limitline_top,'g--');
limitline_bottom = Tmin*ones(1,length(t));
plot(t,limitline_bottom,'g--');
xlabel('t (hours)')
ylabel('Temperature (°C)')
title('Indoor and outdoor temperature')
axis([0 24 0 45])
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', 0:5:45,'LineWidth', 1)
legend('Outdoor temperature','Indoor temperature');
xticks(0:4:24)
xticklabels({'00:00','04:00','08:00','12:00','16:00','20:00','24:00'})


