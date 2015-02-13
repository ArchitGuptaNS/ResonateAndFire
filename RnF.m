function [V, P] = RnF(Iext, w, dt)

b = -15;
% w = 100*pi;
th= 1; 

rs=0.001;

% b = -50;
% w = ErbRateToHz(cfs(1))*2*pi*ones(size(cfs))';
% th= 0.01;

% b = -800;
% w = 400*pi;
% th= 0.29;

[N_ M_] = size(Iext);
V = zeros([N_ M_]);
P = cell([N_ 1]);
z = zeros([N_ 1]);

for t = 1:M_
    z = z + dt*(z.*(b+w*1i) + Iext(:,t));
    V(:, t) = imag(z);
    
    for j=1:N_
        if(V(j,t)>=th)
            P{j} = cat(1, P{j}, t);
            z(j) = 0 + rs*i;
        end
    end
end

disp('Neuron Simulated');

