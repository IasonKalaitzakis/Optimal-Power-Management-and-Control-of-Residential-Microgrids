function cost_equation  = psoClean(numParticles,maxEpochs)


persistent currentEpoch currentPar particlesTable bestFvals neighborIndex pen penalties_table EV_bat_n SoC_initial_ev dt Ev_bat_capacity SoC_target_ev Tmax ...
    Tmin T_in_initial V C rho Uwall Uwindow Fwall Fwindow aw Rsej taf_window SC QEC_max QEC_min COP SOC_EV_penalty_multiplier SOC_EV_target_penalty_multiplier ...
    PPEV_diff_penalty_multiplier offset_variable partial_cost_multiplier SoCmax_ev SoCmin_ev numVariables lb ub Pmin_ev Pmax_ev Pmax_bat Pmin_bat SoC_initial_bat ...
    SoC_target_bat bat_capacity bat_n span CostOfPower cSelf cSocial inertia outdoor_temperature solar_radiation heatcool_sign Zeta Alpha Beta Qin Qsw Qsg U_initial U ...
    QEC_diff_penalty_multiplier Tin_penalty_multiplier T t;

% dataType = 'double';

assignin('base','currentEpoch',currentEpoch);
assignin('base','currentPar',currentPar);
assignin('base','particlesTable',particlesTable);

rng shuffle;
% numParticles = 11;
% maxEpochs = 100;

%%%%%%%%%%MALLON LATHOS, FIX LATER
% minNeighborhoodSize = max(2,floor(numParticles*0.25)); %0.25 =Minimum adaptive neighborhood size, a scalar from 0 to 1. Default is 0.25.
minNeighborhoodSize = 2;
%%%%%%%%%%MALLON LATHOS, FIX LATER

tfInvalidMask = zeros(numParticles,3*numVariables+1);

%Case: First run of function, initialize variables
if isempty(currentEpoch)
    
    numVariables = 49+25;
    cSelf = 1.496;
    cSocial = 1.496;
%     inertia = 0.798;
    inertia = 0.398;        %TODO: make inertia vector for each variable
    
    dt=0.5;
    T=24;
    t = 0:dt:T;
    
    % CostOfPower(1:14)= 0.07897; CostOfPower(15:46)= 0.11058; CostOfPower(47:49)= 0.07897; %€/kwh
    COP_Ther_Var_struct = load('cop_ther_var.mat');
    COP_Ther_Var=COP_Ther_Var_struct.vq1;               %Cost of Power - Summer - Variable
    CostOfPower = COP_Ther_Var;
    
    Pmax_ev = 5;
    Pmin_ev = -5;
    EV_bat_n = 0.85;
    SoC_initial_ev = 19.92;
    dt = 0.5;
    Ev_bat_capacity = 35;
    SoC_target_ev = 20;
    
    Pmax_bat = 5;
    Pmin_bat = -5;
    SoC_initial_bat = 8;
    SoC_target_bat = 8;
    bat_capacity = 15;
    bat_n = 0.85;
    
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
    Alpha = 1/Zeta;
    Beta = (Zeta-(Uwall*Fwall)-(Uwindow*Fwindow))/Zeta;
    Qin = 0.2*ones(1,length(t));
    Qsw = aw * Rsej * Uwall * Fwall * solar_radiation;
    Qsg = 0.001*taf_window * SC * Fwindow * solar_radiation;
%     Qin = zeros(1,length(t));
%     Qsw = zeros(1,length(t));
%     Qsg = zeros(1,length(t));
    
    U = ((Uwall*Fwall+Uwindow*Fwindow)/Zeta)*outdoor_temperature+(Qin+Qsw+Qsg)/Zeta
    U_initial=((Uwall*Fwall+Uwindow*Fwindow)/Zeta)*outdoor_temperature(1);
    
    
    
    SOC_EV_penalty_multiplier = 0.05;
    SOC_EV_target_penalty_multiplier = 1.5;
    PPEV_diff_penalty_multiplier = 0.05;
    
    QEC_diff_penalty_multiplier = 0.2;          % Penalty for large deviations of the QEC
    Tin_penalty_multiplier = 0.6;                 % Penalty for not maintaining the internal temperature between the desired limits
    
    offset_variable = 4; %Offset the sum of power rates to avoid reaching negative values
    partial_cost_multiplier = 8;
    
    SoCmax_ev = Ev_bat_capacity*0.9;
    SoCmin_ev = Ev_bat_capacity*0.1;
    
    lb = [Pmin_ev*ones(1,49) QEC_min*ones(1,25)];
    ub = [Pmax_ev*ones(1,49) QEC_max*ones(1,25)];
    span = ub - lb;
    
    
    currentEpoch = 0;
    currentPar = 0;
    particlesTable = zeros(numParticles,3*numVariables+1); %kathe numVariables plithos timon, gia kathe dianusma (current pos, best pos, velocity)
    particlesTable(:,1:numVariables) = repmat(lb,numParticles,1) + repmat(span,numParticles,1) .* rand(numParticles,numVariables);
    particlesTable(:,numVariables+2:numVariables+2+numVariables-1) = particlesTable(:,1:numVariables);
    % Initialize velocities by randomly sampling over the smaller of options.InitialSwarmSpan or ub-lb. Note that min will be InitialSwarmSpan if either lb or ub is not finite.
    vmax = span;
    particlesTable(:,2*numVariables+2:2*numVariables+2+numVariables-1) = repmat(-vmax,numParticles,1) + ...
        repmat(2*vmax,numParticles,1) .* rand(numParticles,numVariables);
    bestFvals = 0;
    
    
    
    
end

%Case: First iteration, initialize variables
if currentEpoch == 0
    
    currentPar = currentPar + 1;
    
    
    %% ------------- EV Calculations -------------
    PPEV = particlesTable(currentPar,1:49);  % TODO Fix number of variables once you add more power arrays to calculate
    
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
    QEC_not_interp = particlesTable(currentPar,50:74);  % TODO Fix number of variables once you add more power arrays to calculate
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
    
    
    %% ------------- EV Penalty -------------
    SOC_ev_penalty = 0;
    for i=1:length(state_of_charge_ev)
        if state_of_charge_ev(i)>SoCmax_ev
            SOC_ev_penalty = SOC_ev_penalty + (state_of_charge_ev(i) - SoCmax_ev);
        elseif  state_of_charge_ev(i)<SoCmin_ev
            SOC_ev_penalty = SOC_ev_penalty + (SoCmin_ev - state_of_charge_ev(i));
        end
    end
    SOC_ev_target_penalty= abs(state_of_charge_ev(end)-SoC_target_ev);
    
    PPEV_diff_penalty = sum(abs(diff(PPEV)));
    
    
    %% ------------- Heat System Penalty -------------
    Tin_penalty = 0;
    for i=1:length(QEC)
        if indoor_temperature(i)>Tmax(1)
            Tin_penalty = Tin_penalty + (indoor_temperature(i) - Tmax(1));
        elseif  indoor_temperature(i)<Tmin(1)
            Tin_penalty = Tin_penalty + (Tmin(1) - indoor_temperature(i));
        end
    end
    QEC_diff_penalty=sum(abs(diff(QEC)    )   );
    
    
    %% ------------- Cost Calculations -------------
    cost_partial = sum((PPEV+abs(QEC)).*CostOfPower)*dt;    
    cost_equation =(SOC_EV_target_penalty_multiplier * SOC_ev_target_penalty + SOC_EV_penalty_multiplier * SOC_ev_penalty + ...
        PPEV_diff_penalty_multiplier * PPEV_diff_penalty + QEC_diff_penalty_multiplier*QEC_diff_penalty + Tin_penalty_multiplier*Tin_penalty )+ ...
        (partial_cost_multiplier*(offset_variable+cost_partial));
    
    particlesTable(currentPar,numVariables+1) = cost_equation; %fval of previous particle stored as the best value of this particle (since it is the first time that it gets a value)
    
    %Finalize initialization
    if currentPar == numParticles
        currentEpoch = 1;
        currentPar = 0;
        bestFvals = max(particlesTable(:,numVariables+1));
    end
    
elseif currentEpoch <= maxEpochs %main loop
    
    currentPar = currentPar + 1;
    
    if currentPar == 1
        fprintf('Epoch -> %d\n',currentEpoch);
        
        % Generate a random neighborhood for each particle that includes the particle itself
        neighborIndex = zeros(numParticles, minNeighborhoodSize);
        neighborIndex(:, 1) = 1:numParticles; % First neighbor is self
        for i = 1:numParticles
            % Determine random neighbors that exclude the particle itself, which is (numParticles-1) particles
            neighbors = randperm(numParticles-1, minNeighborhoodSize-1);
            % Add 1 to indices that are >= current particle index
            iShift = neighbors >= i;
            neighbors(iShift) = neighbors(iShift) + 1;
            neighborIndex(i,2:end) = neighbors;
        end
        
        % Identify the best neighbor
%         [~, bestRowIndex] = max(particlesTable(neighborIndex(:,2:end)), [], 2); %LATHOS
        % Create the linear index into neighborIndex
%         bestLinearIndex = (bestRowIndex.'-1).*numParticles + (1:numParticles);
        %         bestNeighborIndex = neighborIndex(bestLinearIndex);
        randSelf = rand(numParticles, numVariables);
        randSocial = rand(numParticles, numVariables);
        
        % Update the velocities
        
        temp1 = inertia*particlesTable(:,2*numVariables+2:2*numVariables+2+numVariables-1);
        temp2 = cSelf*randSelf.*(particlesTable(:,numVariables+2:numVariables+2+numVariables-1)-particlesTable(:,1:numVariables));
        %  temp3 =  cSocial*randSocial.*(particlesTable(bestNeighborIndex,numVariables+2:numVariables+2+numVariables-1)-particlesTable(:,1:numVariables))
        
        %         bestLinearIndex
        %         neighborIndex
        %         particlesTable(neighborIndex(:,2:end),numVariables+2:numVariables+2+numVariables-1)
        %         particlesTable(:,1:numVariables)
        
        temp3 =  cSocial*randSocial.*(particlesTable(neighborIndex(:,2:end),numVariables+2:numVariables+2+numVariables-1)-particlesTable(:,1:numVariables));
        newVelocities =  temp1 + ...
            temp2 + ... %%thesi poy eixe kaluteri timi - current thesi
            temp3;
        
        
        %%%         tfValid = all(isfinite(newVelocities), 2);
        %%%         state.Velocities(tfValid,:) = newVelocities(tfValid,:);
        particlesTable(:,2*numVariables+2:2*numVariables+2+numVariables-1) = newVelocities;
        
        % Update the positions
        newPopulation = particlesTable(:,1:numVariables) + particlesTable(:,2*numVariables+2:2*numVariables+2+numVariables-1);
        %%%         tfInvalid = ~isfinite(newPopulation);
        %%%         newPopulation(tfInvalid) = state.Positions(tfInvalid);
        % Enforce bounds, setting the corresponding velocity component to
        % zero if a particle encounters a lower/upper bound
        tfInvalid = newPopulation < lb;
        tfInvalidMask(:,2*numVariables+2:end) = tfInvalid;
        mask = (tfInvalidMask==1);
        
        
        if any(any(tfInvalid))
            %newPopulation(tfInvalid) = lb; %%Replace this command with the following loop due to errors
            for y=1:numVariables
                for i=1:numParticles
                    if tfInvalid(i,y)==1
                        newPopulation(i,y)= lb(y);
                    end
                end
            end
            particlesTable(mask) = 0;
            %             disp('Applied low tfInvalid')
        end
        
        tfInvalid = newPopulation > ub;
        tfInvalidMask(:,2*numVariables+2:end) = tfInvalid;
        mask = (tfInvalidMask==1);
        if any(any(tfInvalid))
            %             newPopulation(tfInvalid) = ub;
            for y=1:numVariables
                for i=1:numParticles
                    if tfInvalid(i,y)==1
                        newPopulation(i,y)= ub(y);
                    end
                end
            end
            particlesTable(mask) = 0;
            %             disp('Applied high tfInvalid')
        end
        particlesTable(:,1:numVariables) = newPopulation;
    end
    
    
    %% ------------- EV Calculations -------------
    PPEV = particlesTable(currentPar,1:49);  % TODO Fix number of variables once you add more power arrays to calculate
    
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
    QEC_not_interp = particlesTable(currentPar,50:74);  % TODO Fix number of variables once you add more power arrays to calculate
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
    
    %% ------------- EV Penalty -------------
    SOC_ev_penalty = 0;
    for i=1:length(state_of_charge_ev)
        if state_of_charge_ev(i)>SoCmax_ev
            SOC_ev_penalty = SOC_ev_penalty + (state_of_charge_ev(i) - SoCmax_ev);
        elseif  state_of_charge_ev(i)<SoCmin_ev
            SOC_ev_penalty = SOC_ev_penalty + (SoCmin_ev - state_of_charge_ev(i));
        end
    end
    SOC_ev_target_penalty= abs(state_of_charge_ev(end)-SoC_target_ev);
    
    PPEV_diff_penalty = sum(abs(diff(PPEV)));
    
    
    %% ------------- Heat System Penalty -------------
    Tin_penalty = 0;
    for i=1:length(QEC)
        if indoor_temperature(i)>Tmax(1)
            Tin_penalty = Tin_penalty + (indoor_temperature(i) - Tmax(1));
        elseif  indoor_temperature(i)<Tmin(1)
            Tin_penalty = Tin_penalty + (Tmin(1) - indoor_temperature(i));
        end
    end
    QEC_diff_penalty=sum(abs(diff(QEC)    )   );
    
    
    %% ------------- Cost Calculations -------------
    cost_partial = sum((PPEV+abs(QEC)).*CostOfPower)*dt;
    %     cost_equation =(1+ SOC_EV_target_penalty_multiplier * SOC_ev_target_penalty + SOC_EV_penalty_multiplier * SOC_ev_penalty + ...
    %         PPEV_diff_penalty_multiplier * PPEV_diff_penalty)*(offset_variable+cost_partial);
    
    cost_equation =(SOC_EV_target_penalty_multiplier * SOC_ev_target_penalty + SOC_EV_penalty_multiplier * SOC_ev_penalty + ...
        PPEV_diff_penalty_multiplier * PPEV_diff_penalty + QEC_diff_penalty_multiplier*QEC_diff_penalty + Tin_penalty_multiplier*Tin_penalty )+ ...
        (partial_cost_multiplier*(offset_variable+cost_partial));
    
    if cost_equation < particlesTable(currentPar,numVariables+1)

        particlesTable(currentPar,numVariables+1) = cost_equation;
        particlesTable(currentPar,numVariables+2:numVariables+2+numVariables-1) = particlesTable(currentPar,1:numVariables);
        
        
        penalties_table(currentPar,1) = SOC_EV_target_penalty_multiplier * SOC_ev_target_penalty;
        penalties_table(currentPar,2) = SOC_EV_penalty_multiplier * SOC_ev_penalty;
        penalties_table(currentPar,3) = PPEV_diff_penalty_multiplier * PPEV_diff_penalty;
        penalties_table(currentPar,4) = QEC_diff_penalty_multiplier * QEC_diff_penalty;
        penalties_table(currentPar,5) = Tin_penalty_multiplier * Tin_penalty;
        penalties_table(currentPar,6) = partial_cost_multiplier*(offset_variable+cost_partial);
        penalties_table(currentPar,7) = particlesTable(currentPar,numVariables+1);

    end
    
    %debugging
    if currentPar == 1
        [minimum,~] = min(particlesTable(:,numVariables+1));
        disp(minimum)
    end
    %debugging
    if currentEpoch == maxEpochs
        pen = array2table(penalties_table);
        pen.Properties.VariableNames(1:7) = {'EV SoC target','EV SoC','EV Diff','QEC Diff','Temperature','offset+cost','Final cost'};
        assignin('base','pen',pen);
    end
    
    
    if currentPar == numParticles
        
        currentEpoch = currentEpoch + 1;
        currentPar = 0;
        return;
        
    elseif currentPar>1 && currentPar<numParticles
        return;
        
    elseif currentPar == 1
        return;
        
    else % never reaches that
        error('Particle Loop Error')
    end
else
    
    
    [cost_equation,i_min_cost] = min(particlesTable(:,numVariables+1));
    assignin('base','i_min_cost',i_min_cost);
    i_min_cost
    
    disp('End of iterations');
    
end
end




