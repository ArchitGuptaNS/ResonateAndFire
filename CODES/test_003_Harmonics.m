close all;
clear all;
clc;

INSTRUMENTS = {'Viola'; 'Violin'; 'Guitar'; 'Ukulele'; 'Bass'; 'Banjo'};
BASE_ADDRS  = '../CASES/';

i           = 1;
j           = 1;

ADDR        = strcat(BASE_ADDRS,INSTRUMENTS{i}, '/', int2str(j), '.wav');
file        = wavread(ADDR);