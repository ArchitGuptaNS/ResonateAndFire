clc;
close all;
clear all;

%%  INITIALIZATION  %%

dt      = 1e-4;
T_MAX   = 1.4;
CF      = 100*pi;
NC      = 8;
SCALE   = 80;

in = SCALE*sin(CF*(dt:dt:T_MAX));
on = in;
in(end-4000:end) = 0;
% in(end-2000:end) = in(end-2000:end)/2;
w  = pow2((-NC:NC)/NC)' * CF;
in = repmat(in, [length(w) 1]);

Wts   = 0.5*ones(2*NC+1);
for ij=1:2*NC+1
   Wts(ij,ij) = 0; 
end

NITERS  = 1;

%%  TESTING BEHAVIOUR FOR NEARBY CENTRAL FREQUENCIES  %%

[V, P]  = RnF(in, w, dt);
figure();
plot(dt:dt:T_MAX, in);
hold on;
plot(dt:dt:T_MAX, SCALE*V(NC+1,:));

figure();
for dm = -4:4
    subplot(3, 3, 5+dm), plot(dt:dt:T_MAX, V(NC+1+dm,:));
    title(int2str(w(NC+1+dm)));
end

%%  INTRODUCING LATERAL CONNECTIONS WITH FIXED WEIGHTS %%

Iext    = in;
b       = -15;
th      = 1; 
rs      =0.001;

Aup = 0.01;
Adn = -0.02;
Tl  = 200;
T   = 150;
ts  = T/4;
I0  = SCALE;

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
                z(j) = 0 + rs*i;
            end
        end
    end
end

%%  PLOTS AND ANALYSIS  %%

tst_001_INPUT = figure();
plot(dt:dt:T_MAX, in);
hold on;
plot(dt:dt:T_MAX, SCALE*V(NC+1,:));

tst_001_Current = figure();
plot(dt:dt:T_MAX, Iext(NC+1,:));
hold on;
plot(dt:dt:T_MAX, in);

tst_001_Weights = figure();
plot(Wts(NC+1,:));

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

figure();
plot(dt:dt:T_MAX, on/SCALE);
hold on;
plot(dt:dt:T_MAX, in/SCALE);
for dP = 1:length(P)
    scatter(P{dP}, ((dP-NC-1)/NC)*ones(size(P{dP})));
    hold on;
end