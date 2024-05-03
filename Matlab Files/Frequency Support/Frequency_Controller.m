function DPdem = Frequency_Controller(f,DPmax,DPmin)

fg = 50; % Nominal frequency of the grid
f1 = 0.984*fg; %49,2
f2 = 0.996*fg; %49,8
f3 = 1.004*fg; %50,2
f4 = 1.016*fg; %50,8

% plot_droop(f1,f2,f3,f4,DPmax,DPmin)

if f < f1
    DPdem = DPmin;
elseif f1 <= f && f <= f2
    DPdem = DPmin - (DPmin/(f2-f1))*(f-f1);
elseif f2 <= f && f <= f3
    DPdem = 0;
elseif f3 <= f && f <= f4
	DPdem = (DPmax/(f4-f3))*(f-f3);
else %f > f4
    DPdem = DPmax;
end

% Pref = P + DPdem; 


