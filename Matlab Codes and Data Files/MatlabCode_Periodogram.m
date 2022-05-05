%---- Code---
clc; clear all; format long;
%load snowD.xlsx
load -ascii SnowDepth_Data_56Days_nsnow.plt
k=1;
for i=1:1344
 t(k)= SnowDepth_Data_56Days_nsnow(i);
 snowD(k)= SnowDepth_Data_56Days_nsnow(1344+i);
 k=k+1;
end
snowDnorm = snowD-mean(snowD);
fs = 24*19;
t = (0:length(snowDnorm) - 1)/fs;
figure(1)
plot(t,snowDnorm)
xlabel('Time (days)')
ylabel('SnowDepth ({}m)')
axis tight


[pxx,f] = periodogram(snowDnorm,[],[],fs);
figure(2)
plot(f,pxx)
ax = gca;
ax.XLim = [0 56];
xlabel('Frequency (Days)')
ylabel('Magnitude (Snow Depth)')
