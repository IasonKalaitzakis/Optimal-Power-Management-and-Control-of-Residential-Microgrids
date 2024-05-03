clear all; close all; clc;

load('kanellos_R1.mat')
p_PEV(:,1)=PPEV;
p_BAT(:,1)=PBAT;
p_QEC(:,1)=QEC_electric;
p_LOAD(:,1)=PLOAD;
socBAT(:,1) = state_of_charge_bat;
socPEV(:,1) = state_of_charge_ev;
indoor_temperature_table(:,1)=indoor_temperature;
outdoor_temperature(:,1)=outdoor_temperature_temp;
home_load(:,1)=home_load_temp;
PV_Production(:,1)=PV_Production_temp;
cost(1)=final_cost;

load('kanellos_R2.mat')
p_PEV(:,2)=PPEV;
p_QEC(:,2)=QEC_electric;
p_LOAD(:,2)=PLOAD;
socPEV(:,2) = state_of_charge_ev;
indoor_temperature_table(:,2)=indoor_temperature;
outdoor_temperature(:,2)=outdoor_temperature_temp;
home_load(:,2)=home_load_temp;
cost(2)=final_cost;

load('kanellos_R3.mat')
p_BAT(:,3)=PBAT;
p_LOAD(:,3)=PLOAD;
socBAT(:,3) = state_of_charge_bat;
home_load(:,3)=home_load_temp;
PV_Production(:,3)=PV_Production_temp;
cost(3)=final_cost;

load('kanellos_R4.mat')
p_PEV(:,4)=PPEV;
p_QEC(:,4)=QEC_electric;
p_LOAD(:,4)=PLOAD;
socPEV(:,4) = state_of_charge_ev;
indoor_temperature_table(:,4)=indoor_temperature;
outdoor_temperature(:,4)=outdoor_temperature_temp;
home_load(:,4)=home_load_temp;
cost(4)=final_cost;

p_GRID(:,1) = p_PEV(:,1) + p_BAT(:,1) + p_QEC(:,1) + p_LOAD(:,1) - PV_Production(:,3);
p_GRID(:,2) = p_PEV(:,2) + p_QEC(:,2) + p_LOAD(:,2);
for i=1:49
    if p_GRID(i,2) > 10
        p_GRID(i,2) = 10;
    end
end
p_GRID(:,3) = p_BAT(:,3) + p_LOAD(:,3) - PV_Production(:,3);
p_GRID(:,4) = p_PEV(:,4) + p_QEC(:,4) + p_LOAD(:,4);

p_total_GRID = sum(p_GRID,2);


figure
hold on;
plot(t,p_GRID(:,1),'b','LineWidth', 1.5);
plot(t,p_GRID(:,2),'r','LineWidth', 1.5);
plot(t,p_GRID(:,3),'g','LineWidth', 1.5);
plot(t,p_GRID(:,4),'m','LineWidth', 1.5);
xBox = [0, 24, 24, 0];
yBox = [0, 0, 10, 10];
patch(xBox, yBox, 'black', 'FaceColor', 'green', 'FaceAlpha', 0.1);
xBox = [0 24 24 0];
yBox = [0 0 -10 -10];
patch(xBox, yBox, 'black', 'FaceColor', 'red', 'FaceAlpha', 0.1);
axis([0 24 -10 10])
title('Grid Power exchange of each residence')
xlabel('t (hours)')
ylabel('P_g_r_i_d (kW)')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', -10:2.5:10,'LineWidth', 1)
set(gca, 'FontName', 'Helvetica')
xticks(0:4:24)
xticklabels({'00:00','04:00','08:00','12:00','16:00','20:00','24:00'})
legend('Residence 1','Residence 2','Residence 3','Residence 4');


figure
hold on;
plot(t,p_total_GRID,'b','LineWidth', 1.5);
xBox = [0, 24, 24, 0];
yBox = [0, 0, 30, 30];
patch(xBox, yBox, 'black', 'FaceColor', 'green', 'FaceAlpha', 0.1);
xBox = [0 24 24 0];
yBox = [0 0 -30 -30];
patch(xBox, yBox, 'black', 'FaceColor', 'red', 'FaceAlpha', 0.1);
axis([0 24 -30 30])
title('Total Microgrid Power Exchange')
xlabel('t (hours)')
ylabel('P_g_r_i_d (kW)')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', -30:5:30,'LineWidth', 1)
set(gca, 'FontName', 'Helvetica')
xticks(0:4:24)
xticklabels({'00:00','04:00','08:00','12:00','16:00','20:00','24:00'})




