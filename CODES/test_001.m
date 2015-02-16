clc;
close all;
clear all;

%%  INITIALIZATION  %%

dt      = 1e-4;
T_MAX   = 1.4;
CF      = 100*pi;
NC      = 12;
SCALE   = 80;
NOISE_LEVEL = 1;

on = SCALE*sin(CF*(dt:dt:T_MAX));
in = on;
Iext= [];
dec= pow2((0:-1:-4200)/1000);
w  = pow2((-NC:NC)/(1.2*NC))' * CF;

for ij=-NC:NC
    Iext = cat(1, Iext, SCALE*sin(CF*(dt:dt:T_MAX)+abs(ij/NC)));
%     Iext(NC+1+ij,end-4200:end) = Iext(NC+1+ij,end-4200:end) .* dec;
    Iext(NC+1+ij,end-4200:end) = 0;
end

Iext = Iext + NOISE_LEVEL*randn(size(Iext));
in = Iext(NC+1,:);

Wts   = 0.3*ones(2*NC+1);

for ij=1:2*NC+1
   Wts(ij,ij) = 0; 
end

NITERS  = 1;

%%  TESTING BEHAVIOUR FOR NEARBY CENTRAL FREQUENCIES  %%

[V, P]  = RnF(Iext, w, dt);

tst_001_INPUT = figure();
subplot(2,1,1), plot(dt:dt:T_MAX, in, 'color', 'red');
hold on;
plot(dt:dt:T_MAX, SCALE*V(NC+1,:));
legend('Input Current', 'Output Voltage');
title('WITHOUT LATERAL CONNECTIONS');
xlim([.4 1.2]);

test_001_hist = figure();
subplot(2,1,1);
plot(dt:dt:T_MAX, in/SCALE);
title('WITHOUT LATERAL CONNECTIONS');
hold on;
for dP = 1:length(P)
    scatter(P{dP}*dt, ((dP-NC-1)/NC)*ones(size(P{dP})));
    hold on;
end

figure();
for dm = -4:4
    subplot(3, 3, 5+dm), plot(dt:dt:T_MAX, V(NC+1+dm,:));
    title(int2str(w(NC+1+dm)));
end

%%  INTRODUCING LATERAL CONNECTIONS WITH FIXED WEIGHTS %%

b       = -15;
th      = 1; 
rs      =0.01;

Aup = 0.2;
Adn = -0.02;
Tl  = 200;
T   = 150;
ts  = T/4;
I0  = 1.28*SCALE;

for dIter = 1:NITERS
    
    [N_ M_] = size(Iext);
    V = zeros([N_ M_]);
    P = cell([N_ 1]);
    z = zeros([N_ 1]);

    for t = 1:M_
        z = z + dt*(z.*(b+w*1i) + Iext(:,t));
        V(:, t) = imag(z);

        for j=1:N_
            if(V(j,t)>=th)
                P{j} = cat(1, P{j}, t*dt);
                tm   = (t+1 : M_) - t;
                Iext(:,t+1:end) = Iext(:,t+1:end) + I0 * Wts(:,j) * (exp(-tm/T)-exp(-tm/ts));
                for jk = 1:N_
                    if(~isempty(P{jk}))
                        Wts(jk,j) = Wts(jk,j)*(1+Aup*(exp((P{jk}(end) - t)/Tl)));
                        Wts(j,jk) = Wts(j,jk)*(1+Adn*(exp((P{jk}(end) - t)/Tl)));
                    end
                end
                z(j) = rs*z(j);
            end
        end
    end
end

%%  PLOTS AND ANALYSIS  %%
figure(tst_001_INPUT);
subplot(2,1,2), plot(dt:dt:T_MAX, in, 'color', 'red');
hold on;
plot(dt:dt:T_MAX, SCALE*V(NC+1,:));
legend('Input Current', 'Output Voltage');
title('WITH LATERAL CONNECTIONS');
xlim([0.4 1.2]);

tst_001_Current = figure();
subplot(2,1,1);
plot(dt:dt:T_MAX, Iext(NC+1,:), 'color', [0.7 0.7 0.7]);
hold on;
plot(dt:dt:T_MAX, in);
legend('Current POST STDP', 'Original Current');

subplot(2,1,2);
plot(dt:dt:T_MAX, Iext(NC+1,:)-in); 
title('Difference between the two');

tst_001_Weights = figure();
plot(Wts(NC+1,:));
title('Synaptic Weights');

tst_001_OUTPUT = figure();
SHOW = 5;
SRT  = ceil((SHOW*SHOW - 1)/2);
for dm = -SRT:SRT
    if(abs(dm) <= NC)
        subplot(SHOW, SHOW, SRT+1+dm), plot(dt:dt:T_MAX, V(NC+1+dm,:));
        hold on;
        plot(dt:dt:T_MAX, Iext(NC+1+dm,:)/I0, 'color', 'red');
        title(int2str(w(NC+1+dm)));
        xlim([0.6 1.2]);
    end
end

figure(test_001_hist);
title('RASTER PLOT FOR SPIKE TIMINGS');
subplot(2,1,2);
plot(dt:dt:T_MAX, in/SCALE);
title('WITH LATERAL CONNECTIONS');
hold on;
for dP = 1:length(P)
    scatter(P{dP}, ((dP-NC-1)/NC)*ones(size(P{dP})));
    hold on;
end