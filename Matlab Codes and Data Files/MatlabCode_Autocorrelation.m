%---- Code---
clc; clear all; format long;
%load snowD.xlsx
load -ascii SnowDepth_Data_56Days_nsnow.plt
k=1;
for i=1:1344
 t(k)= (i);
 snowD(k)= SnowDepth_Data_56Days_nsnow(1344+i);
 k=k+1;
end
snowDnorm = snowD-mean(snowD);
fs = 24;
t = (0:length(snowDnorm) - 1)/fs;
figure(1)
plot(t,snowDnorm)
xlabel('Time (days)')
ylabel('SnowDepth ({}m)')
axis tight
%%%%%
[autocor,lags] = xcorr(snowDnorm,100*fs,'coeff');
figure(2)
plot(lags/fs,autocor)
xlabel('Lag (days)')
ylabel('Autocorrelation')
axis([-30 30 -0.4 1])
%%%%
[pksh,lcsh] = findpeaks(autocor);
short = mean(diff(lcsh))/fs
%%%%%%
[pklg,lclg] = findpeaks(autocor, ...
    'MinPeakDistance',ceil(short)*fs,'MinPeakheight',0.04);
long = mean(diff(lclg))/fs
%%%%%%%
figure(3)
hold on
pks = plot(lags(lcsh)/fs,pksh,'or');%, ...
   % lags(lclg)/fs,pklg+0.05,'vk');
hold off
%legend(pks,[repmat('Period: ',[2 1]) num2str([short;long],0)])
xlabel('Lag (days)')
ylabel('Autocorrelation')
axis([-30 30 -0.4 1])



