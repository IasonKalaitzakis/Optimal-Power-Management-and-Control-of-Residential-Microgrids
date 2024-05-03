function pso_out=prepare_PSO( CostOfPower_t,  CostOfPV_t,  PV_Production_t,  home_load_t,  SoCmax_ev_t,  SoCmin_ev_t,  SoC_initial_ev_t,  SoC_commute_t,  Tmax_t,  Tmin_t,  T_in_initial_t,...
                 Zeta_t,  SoC_target_ev_t,  dt_t,  outdoor_temperature_t,  Uwall_t,  Uwindow_t,  Fwall_t,  Fwindow_t,  Qin_t,  Qsw_t,  Qsg_t,  dim1_t,  dim2_t,  dim3_t,  dim4_t,...
                 heatcool_sign_t,  t_leave_t,  t_return_t,  SoCmax_bat_t,  SoCmin_bat_t,  SoC_initial_bat_t,  SoC_target_bat_t,  EV_bat_n_t,  bat_n_t,...
                functname, dims, mv, varrange, minmax, psoparams,plotfcn,OptimisationStartTime_t,EV_flag_t,Battery_flag_t,HeatSystem_flag_t,Loadshift_flag_t,COP_t,commute_enable_t,grid_outage_t,SellingCostOfPower_t)
            
%The purpose of this code is to pass the globals into PSOfunction for the
%GUI use, since app design does not support globals
            
global CostOfPower CostOfPV PV_Production home_load SoCmax_ev SoCmin_ev SoC_initial_ev ...
SoC_commute Tmax Tmin T_in_initial Zeta SoC_target_ev dt outdoor_temperature Uwall Uwindow Fwall ...
Fwindow Qin Qsw Qsg dim1 dim2 dim3 dim4 heatcool_sign t_leave t_return SoCmax_bat SoCmin_bat ...
SoC_initial_bat SoC_target_bat EV_bat_n bat_n OptimisationStartTime EV_flag Battery_flag ...
HeatSystem_flag Loadshift_flag COP commute_enable grid_outage SellingCostOfPower

CostOfPower=CostOfPower_t;
CostOfPV=CostOfPV_t;
PV_Production=PV_Production_t;
home_load=home_load_t;
SoCmax_ev=SoCmax_ev_t;
SoCmin_ev=SoCmin_ev_t;
SoC_initial_ev=SoC_initial_ev_t;
SoC_commute=SoC_commute_t;
Tmax=Tmax_t;
Tmin=Tmin_t;
T_in_initial=T_in_initial_t;
Zeta=Zeta_t;
SoC_target_ev=SoC_target_ev_t;
dt=dt_t;
outdoor_temperature=outdoor_temperature_t;
Uwall=Uwall_t;
Uwindow=Uwindow_t;
Fwall=Fwall_t;
Fwindow=Fwindow_t;
Qin=Qin_t;
Qsw=Qsw_t;
Qsg=Qsg_t;
dim1=dim1_t;
dim2=dim2_t;
dim3=dim3_t;
dim4=dim4_t;
heatcool_sign=heatcool_sign_t;
t_leave=t_leave_t;
t_return=t_return_t;
SoCmax_bat=SoCmax_bat_t;
SoCmin_bat=SoCmin_bat_t;
SoC_initial_bat=SoC_initial_bat_t;
SoC_target_bat=SoC_target_bat_t;
EV_bat_n=EV_bat_n_t;
bat_n=bat_n_t;
OptimisationStartTime = OptimisationStartTime_t;
EV_flag = EV_flag_t;
Battery_flag = Battery_flag_t;
HeatSystem_flag = HeatSystem_flag_t;
Loadshift_flag = Loadshift_flag_t;
COP = COP_t;
commute_enable = commute_enable_t;
grid_outage = grid_outage_t;
SellingCostOfPower = SellingCostOfPower_t;

[pso_out,~,~]=pso_Trelea_vectorized(functname, dims, mv, varrange, minmax, psoparams, plotfcn);

