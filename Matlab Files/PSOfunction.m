  function cost =PSOfunction(x)
s=size(x);

% [CostOfPower,home_load,SoCmax_ev,SoCmin_ev,SoC_initial_ev,...
%     SoC_commute,Tmax,Tmin,T_in_initial,Zeta,SoC_target_ev,dt,outdoor_temperature,Uwall,Uwindow,Fwall,...
%     Fwindow,Qin,Qsw,Qsg,dim1,dim2,dim3,dim4,heatcool_sign,t_leave,t_return,SoCmax_bat,SoCmin_bat,...
%     SoC_initial_bat,SoC_target_bat,EV_bat_n,bat_n,OptimisationStartTime,EV_flag,Battery_flag,...,
%     HeatSystem_flag,Loadshift_flag,COP,commute_enable,ps,epoch,solar_panels_state,CostOfPV,solar_panel_power,...
%     Pmax_bat,Pmin_bat,Pmax_ev,Pmin_ev,QEC_max,QEC_min,Load_shifting,Alpha,Beta,U_initial,U,scenario]=preparePIL();

global CostOfPower CostOfPV PV_Production home_load SoCmax_ev SoCmin_ev SoC_initial_ev ...
SoC_commute Tmax Tmin T_in_initial Zeta SoC_target_ev dt outdoor_temperature Uwall Uwindow Fwall ...
Fwindow Qin Qsw Qsg dim1 dim2 dim3 dim4 heatcool_sign t_leave t_return SoCmax_bat SoCmin_bat ...
SoC_initial_bat SoC_target_bat EV_bat_n bat_n OptimisationStartTime EV_flag Battery_flag HeatSystem_flag Loadshift_flag COP commute_enable grid_outage SellingCostOfPower




T_end = 24-((OptimisationStartTime-1)/2);

%% ------------- Penalty Multipliers -------------
PPEV_diff_penalty_multiplier = 0.1;        % Penalty for large deviations of the PPEV
PBAT_diff_penalty_multiplier = 0.1;         % Penalty for large deviations of the PBAT

SOC_EV_penalty_multiplier = 0.5;            % Penalty for crossing the imposed SoC limits of the EV battery
SOC_EV_target_penalty_multiplier = 5;     % Penalty for not reaching the SoC target

QEC_diff_penalty_multiplier = 0.02;          % Penalty for large deviations of the QEC
Tin_penalty_multiplier = 3;                 % Penalty for not maintaining the internal temperature between the desired limits

SOC_BAT_penalty_multiplier = 5;           % Penalty for crossing the imposed SoC limits of the battery & for not reaching the SoC target

LOAD_equals_penalty_multiplier = 3;         % Penalty for the initial and load shifted loads not having equal sums
grid_outage_multiplier = 1;

% partial_cost_multiplier = 18.5;
partial_cost_multiplier = 18;

offset_variable =0; %Offset the sum of power rates to avoid reaching negative values



for ii=1:s(1)
    
    SOC_ev_penalty = 0; PPEV_diff_penalty=0; SOC_bat_penalty = 0; PBAT_diff_penalty=0;
    Tin_penalty = 0; QEC_diff_penalty=0; Load_equals_penalty=0; grid_outage_penalty = 0;
    %% ------------- EV Calculations -------------
    if EV_flag=="On"
        PPEV_not_interp = x(ii,1:dim1);
        PPEV = interp1([0:T_end/(dim1-1):T_end],PPEV_not_interp,[0:dt:T_end],'linear');
        
        if commute_enable=="On"
            PPEV(2*t_leave:2*t_return) = 0;
        end
        
        state_of_charge_ev = zeros(1,length(PPEV));
        PPEV_actual = zeros(1,length(PPEV))';
        for z=1:length(PPEV_actual)
            if PPEV(z)>=0
                PPEV_actual(z)=PPEV(z)*EV_bat_n;
            else
                PPEV_actual(z)=PPEV(z)/EV_bat_n;
            end
        end
        
        
        state_of_charge_ev(1)=SoC_initial_ev;
        state_of_charge_ev(2:end)=SoC_initial_ev+cumsum(PPEV_actual(1:end-1))*dt;
        
        if commute_enable=="On"
            for i= 1:2*t_return-2*t_leave
                state_of_charge_ev(t_leave+i)=state_of_charge_ev(2*t_leave+i)-i/((2*t_return-2*t_leave)/SoC_commute);
            end
            state_of_charge_ev(2*t_return+1:end)=state_of_charge_ev(2*t_return+1:end)-SoC_commute;
        end
        
    elseif EV_flag=="Off"
        PPEV=zeros(1,length(0:dt:T_end));
        state_of_charge_ev=zeros(1,length(0:dt:T_end));
    end

    %% ------------- Heat System Calculations -------------
    
    if HeatSystem_flag=="On"
        QEC_not_interp = x(ii,dim1+1:dim1+dim2);
        QEC = interp1([0:T_end/(dim2-1):T_end],QEC_not_interp,[0:dt:T_end],'linear');
        
        indoor_temperature=zeros(1,length(QEC));
        indoor_temperature(1) =T_in_initial;
        
        for i=2:length(QEC)
            
            indoor_temperature(i) = indoor_temperature(i-1) + heatcool_sign*QEC(i-1)/Zeta + ...
                (Uwall*Fwall*(outdoor_temperature(i-1)-indoor_temperature(i-1))+ ...
                Uwindow*Fwindow*(outdoor_temperature(i-1)-indoor_temperature(i-1))+ ...
                Qin(i-1)+Qsw(i-1)+Qsg(i-1))/Zeta;
        end
        
        QEC_electric = QEC/COP;
        
    elseif HeatSystem_flag=="Off"
        QEC = zeros(1,length(0:dt:T_end));
        QEC_electric=zeros(1,length(0:dt:T_end));
        indoor_temperature=[];
    end
    
    %% ------------- Load Shifting Calculations -------------
    if Loadshift_flag=="On"
        PLOAD_shifting_not_interp = x(ii,dim1+dim2+1:dim1+dim2+dim3);
        PLOAD_shifting = interp1([0:T_end/(dim3-1):T_end],PLOAD_shifting_not_interp,[0:dt:T_end],'linear')';
        PLOAD = (PLOAD_shifting.*home_load')';
    elseif Loadshift_flag=="Off"
        PLOAD = home_load;
    end
    
    %% ------------- Battery Calculations ------------- 
    if Battery_flag=="On"
        PBAT_not_interp = x(ii,dim1+dim2+dim3+1:dim1+dim2+dim3+dim4);
        PBAT = interp1([0:T_end/(dim4-1):T_end],PBAT_not_interp,[0:dt:T_end],'linear');
        
        state_of_charge_bat= zeros(1,length(PBAT));
        PBAT_actual = zeros(1,length(PBAT))';
        for z=1:length(PBAT_actual)
            if PBAT(z)>=0
                PBAT_actual(z)=PBAT(z)*bat_n;
            else
                PBAT_actual(z)=PBAT(z)/bat_n;
            end
        end
        state_of_charge_bat(:)=SoC_initial_bat+cumsum(PBAT_actual)*dt;

    elseif Battery_flag=="Off"
        PBAT = zeros(1,length(0:dt:T_end));
        state_of_charge_bat = zeros(1,length(0:dt:T_end));
    end
    
    %% ------------- EV Penalty -------------
    SOC_ev_target_penalty=0;
    if EV_flag=="On"
        for i=1:length(PPEV)
            if state_of_charge_ev(i)>SoCmax_ev
                SOC_ev_penalty = SOC_ev_penalty + (state_of_charge_ev(i) - SoCmax_ev);
            elseif  state_of_charge_ev(i)<SoCmin_ev
                SOC_ev_penalty = SOC_ev_penalty + (SoCmin_ev - state_of_charge_ev(i));
            end
        end
%         SOC_ev_penalty = SOC_ev_penalty + abs(state_of_charge_ev(end)-SoC_target_ev);
        SOC_ev_target_penalty= abs(state_of_charge_ev(end)-SoC_target_ev);
        
        if commute_enable=="On"
            PPEV_diff_penalty = sum(abs(diff(PPEV(1:t_leave-1)))) + sum(abs(diff(PPEV(t_return+1:end))));
        else
            PPEV_diff_penalty = sum(abs(diff(PPEV)));
        end
    end
    
    %% ------------- Battery Penalty -------------
    if Battery_flag=="On"
        for i=1:length(PBAT)
            if state_of_charge_bat(i)>SoCmax_bat
                SOC_bat_penalty = SOC_bat_penalty + (state_of_charge_bat(i) - SoCmax_bat);
            elseif  state_of_charge_bat(i)<SoCmin_bat
                SOC_bat_penalty = SOC_bat_penalty + (SoCmin_bat - state_of_charge_bat(i));
            end
        end
        SOC_bat_penalty = SOC_bat_penalty + abs(state_of_charge_bat(end)-SoC_target_bat);
        PBAT_diff_penalty=sum(abs(diff(PBAT)    )   );
    end

    %% ------------- Heat System Penalty -------------
    if HeatSystem_flag=="On"
        for i=1:length(QEC)
            if indoor_temperature(i)>Tmax(1)
                Tin_penalty = Tin_penalty + (indoor_temperature(i) - Tmax(1));
            elseif  indoor_temperature(i)<Tmin(1)
                Tin_penalty = Tin_penalty + (Tmin(1) - indoor_temperature(i));
            end
        end
        QEC_diff_penalty=sum(abs(diff(QEC)    )   );
    end
    
    %% ------------- Load Shifting Penalty -------------
    if Loadshift_flag=="On"
        Load_equals_penalty = abs(sum(PLOAD)-sum(home_load));
    end
    %% ------------- Grid outage Penalty -------------
    grid_outage_penalty = 0;
    if grid_outage==1
        for i = 31:38
            grid_outage_penalty = grid_outage_penalty + abs(sum(PPEV(i)+PBAT(i)+abs(QEC_electric(i))+PLOAD(i)-PV_Production(i)));
        end
    else
        grid_outage_penalty = 0;
    end
    

    
    %% ------------- Final Cost -------------

%     cost_partial = sum((PPEV+PBAT+abs(QEC_electric)+PLOAD).*CostOfPower)*dt
    power_sum = PPEV+PBAT+abs(QEC_electric)+PLOAD-PV_Production;
    cost_partial = 0; 
    for p = 1 : 49
        if power_sum(p)>=0
            cost_partial = cost_partial + (power_sum(p) * CostOfPower(p));
        else
            cost_partial = cost_partial + (power_sum(p) * SellingCostOfPower(p));
        end
    end

    
%     cost_partial = sum (power_sum .* CostOfPower);

    cost(ii,1) =(PPEV_diff_penalty_multiplier*PPEV_diff_penalty+QEC_diff_penalty_multiplier*QEC_diff_penalty + ... 
        PBAT_diff_penalty_multiplier*PBAT_diff_penalty+ SOC_EV_target_penalty_multiplier*SOC_ev_target_penalty + ... 
        SOC_EV_penalty_multiplier*SOC_ev_penalty+SOC_BAT_penalty_multiplier*SOC_bat_penalty + ... 
        Tin_penalty_multiplier*Tin_penalty + LOAD_equals_penalty_multiplier*Load_equals_penalty + grid_outage_multiplier * grid_outage_penalty) + ... 
        (partial_cost_multiplier*(offset_variable+cost_partial));

%     
%     cost_equation =(SOC_EV_target_penalty_multiplier * SOC_ev_target_penalty + SOC_EV_penalty_multiplier * SOC_ev_penalty + ...
%         PPEV_diff_penalty_multiplier * PPEV_diff_penalty + QEC_diff_penalty_multiplier*QEC_diff_penalty + Tin_penalty_multiplier*Tin_penalty )+ ...
%         (partial_cost_multiplier*(offset_variable+cost_partial));
%     
%     
    
end

% test1 = SOC_EV_penalty_multiplier*SOC_ev_penalty;
% test2 = Tin_penalty_multiplier*Tin_penalty;
% test3 = LOAD_equals_penalty_multiplier*Load_equals_penalty;
% test4 = SOC_BAT_penalty_multiplier*SOC_bat_penalty;
% test5 = SOC_EV_target_penalty_multiplier*SOC_ev_target_penalty;
% test6 = PPEV_diff_penalty_multiplier*PPEV_diff_penalty;
% test7 = QEC_diff_penalty_multiplier*QEC_diff_penalty;
% test8 = PBAT_diff_penalty_multiplier*PBAT_diff_penalty;
% test9 = partial_cost_multiplier*(offset_variable+cost_partial);
% test10 = grid_outage_multiplier * grid_outage_penalty;
% disp('  ppev_diff  qec_diff pbat_diff   soc_ev   indtemp  loadequals soc_bat EV_target partial_cost gridoutage' )
% penalties = [test6 test7 test8 test1 test2 test3 test4 test5 test9 test10];
% disp(penalties)

