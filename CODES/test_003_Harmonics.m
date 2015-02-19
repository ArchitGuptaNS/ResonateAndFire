close all;
clear all;
clc;

%%  CONSTANT DECLARATIONS

AIFF_SCALE_FACTOR   = .001;

%%  EXTRACTING MUSICAL DATA AND RESAMPLING

base_addr   = '../CASES/';
instruments = dir(base_addr);
instruments = instruments(3:end); 

i           = 1;
j           = 1;

pieces 	    = dir(strcat(base_addr,instruments(i).name)); 
pieces      = pieces(3:end);
 	
addr        = strcat(base_addr,'/',instruments(i).name,'/',pieces(j).name);
file        = double(aiffread(addr))*AIFF_SCALE_FACTOR;

FS      = 441;
FC      = 100;

S       = resample(file, FC, FS);

