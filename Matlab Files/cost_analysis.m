close all


idx = 8;



% plot(final_array(idx,2,:),final_array(idx,3,:))
t(1:49) = final_array(idx,2,:);
pev(1:49) = final_array(idx,3,:);
socev(1:49) = final_array(idx,4,:);
pbat(1:49) = final_array(idx,5,:);
socbat(1:49) = final_array(idx,6,:);

figure
plot(t,pev)

figure
plot(t,socev)

figure
plot(t,pbat)
figure
plot(t,socbat)