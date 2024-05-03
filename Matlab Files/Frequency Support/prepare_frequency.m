clear all
clc
close all

f_struct = load('DF.mat')
f = f_struct.Df
fault_duration = 10*60; %Seconds
dt = 1/3600;   % hours
t = [0:dt*60:fault_duration/60];

t_old_f = linspace(0,fault_duration/60,100);

f = interp1(t_old_f,f(1:100),t,'linear');
f = f*65; f = f + 50;



plot(t,f)



%
% a = out.simout.data;
% a = 50*a;
% a(1) = 0;
% a(2) = -1;
% a(3) = -2.5;
% a(4) = -4;
% a(5) = -5.5;
% a(6) = -7;
% a(7) = -8.5;
% a(8) = -10;
%
%
% fault_duration = 20*60; %Seconds
% dt = 1/3600;   % hours
% t = [0:dt*60:fault_duration/60];
% t_old = linspace(0,20,size(a,1));
%
% f = interp1(t_old,a,t,'linear');
%
% f = movmean(f,40);
% f = f+1400;f = f*0.01; f = f + 36; %Negative frequency
%
% % f(1:60) = movmean(f(1:60),40)
% % f(61:99) = movmean(f(61:99),40)
% %
% % f(100:400) = movmean(f(100:400),50)
% figure
% plot(f)
% % plot(t,f)
% save('f.mat','f')
