clc;
close all;
clear all;

%%  INITIALIZATION  %%

dt      = 1e-4;
T_MAX   = 1.4;
CF      = 100*pi*pow2(-2:2)';
NC      = 10;
scale   = [20 36 11 65 12];

in = sin(CF*(dt:dt:T_MAX));
on = scale * in;
on(end-4000:end) = 0;

Wts   = 0.5*ones(2*NC+1);
for ij=1:2*NC+1
   Wts(ij,ij) = 0; 
end

NITERS  = 1;

%%  TESTING BEHAVIOUR FOR DIFFERENT HARMONIC CENTRAL FREQUENCIES  %%


%%  PLOTS AND ANALYSIS

figure(), plot(on);