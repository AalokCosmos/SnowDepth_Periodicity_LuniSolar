% The Hurst exponent
%--------------------------------------------------------------------------
% The first 20 lines of code are a small test driver.
% You can delete or comment out this part when you are done validating the 
% function to your satisfaction.
%
% Bill Davidson, quellen@yahoo.com
% 13 Nov 2005
 clc; clear all; format long;
 load -ascii SnowDepth_Data_700Days_NsnowN.plt
k=1;
for i=1:16800
 t(k)= SnowDepth_Data_700Days_NsnowN(i);
 A(k)= SnowDepth_Data_700Days_NsnowN(16800+i);
 k=k+1;
end
%[]=hurst_exponent()
disp('testing Hurst calculation');

n=16800;
data=A ;%rand(1,n)
plot(data);

hurst=estimate_hurst_exponent(data);

[s,err]=sprintf('Hurst exponent = %.2f',hurst);disp(s);

%[hurst] = estimate_hurst_exponent(data0) ;
%--------------------------------------------------------------------------
% This function does dispersional analysis on a data series, then does a 
% Matlab polyfit to a log-log plot to estimate the Hurst exponent of the 
% series.
%
% This algorithm is far faster than a full-blown implementation of Hurst's
% algorithm.  I got the idea from a 2000 PhD dissertation by Hendrik J 
% Blok, and I make no guarantees whatsoever about the rigor of this approach
% or the accuracy of results.  Use it at your own risk.
%
% Bill Davidson
% 21 Oct 2003

