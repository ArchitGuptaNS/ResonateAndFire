close all;
clear all;
clc;

%%  CONSTANT DECLARATIONS

AIFF_SCALE_FACTOR   = 2;

%%  EXTRACTING MUSICAL DATA AND RESAMPLING

base_addr   = '../CASES/';
instruments = dir(base_addr);
instruments = instruments(3:end); 

i           = 3;
j           = 2;

pieces 	    = dir(strcat(base_addr,instruments(i).name)); 
pieces      = pieces(3:end);
 	
addr        = strcat(base_addr,'/',instruments(i).name,'/',pieces(j).name);
file        = double(aiffread(addr))*AIFF_SCALE_FACTOR;

FS      = 441;
FC      = 100;
dt      = 1e-2/FC;

S       = resample(file(:,1), FC, FS)';
% f       = 27.5;             %   A0
% f       = 30.87;            %   B0
% f       = 29.14;            %   Bb0   

H       = 2*pi*pow2(0:4)';
Imat    = repmat(S, size(H));
[V, P]  = RnF(Imat,f*H,dt);

%% PLOTS AND ANALYSIS

figure();
for i=1:5
    subplot(3,2,i);
    plot(V(i,:), 'color', [0.7 0.7 0.7]);
    hold on;
    plot(S/max(S));
    xlim([41000 42000]);
    legend('Neural Response', 'Audio Input');
end