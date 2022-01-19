% Econ 282B Homework 1
% Ekaterina Gurkova, Christopher Saw
% Winter 2022

clear;
%% Parameters
global R beta delta rho tau q dz fe f fx fxs N theta z_min z_max z_n z_x
R = 1.05;
beta = 1/R;
delta = 0.005;
rho = 5;
tau = 1.231;
q = 0.5;
dz = 0.25;
fe = 1;
f = 0.1;
fx = 1.4;
fxs = 10;

% Distribution of productivity
N = 50;
z_min = 0.25;
z_max = z_min+dz*N;
theta = 5;


% Grid for productivity
z = linspace(z_min-dz, z_max+dz, N+2);

% Initial guess for market demand index
%Pi_d_min = fx/(exp(z_min)*tau^(1-rho));
%Pi_d_max = (fx+fxs)/(exp(z_max)*tau^(1-rho));
Pi_d_min = 0;
Pi_d_max = 5;
Pi_d = (Pi_d_min+Pi_d_max)/2;


%% Iterating over Pi_d
max_iter = 100000;
iter = 1;
tol = 0.0001;
err = 1;
abs1 = 1;
abs2 = 1;


for i=1:N+2
    G(i)=(theta*z_min)^theta / z(i)^(theta+1);
end

Vn = linspace(0,0,N);
Vx = linspace(0,0,N);
Vn_new = linspace(0,0,N);
Vx_new = linspace(0,0,N);
while (abs(err) >= tol) && (iter <= max_iter)
    disp(iter);
    % Value function iteration
    while ((abs1 >= tol) || (abs2 >= tol)) && (iter <= max_iter)
        % Solving for productivity cutoffs for exporting/non-exporting
        z_n = log(fx) - log(Pi_d) - (1-rho)*log(tau);
        z_x = log(fx+fxs) - log(Pi_d) - (1-rho)*log(tau);

        % Distribution of productivity
        %pn = (dz/(z_max-z_min)) / (1-(z_n-z_min)/(z_max-z_min));
        %px = (dz/(z_max-z_min)) / ((z_x-z_min)/(z_max-z_min));
        pn = (-(z_min/(z_n+dz))^theta+((z_min/(z_n))^theta))...
            /((z_min/z_n)^theta);
        px = (-(z_min/(z_x))^theta+((z_min/(z_x - dz))^theta))...
            /(1-(z_min/z_x)^theta);


        for i = 2:N-1
            Vx_new(i) = (1+tau^(1-rho))*Pi_d*exp(z(i)) - f - fx -...
                fxs + beta*(1-delta)*(q*(Vx(i+1) + fxs) + ...
                (1-q)*((1-pn)*(Vx(i-1)+fxs)+pn*Vn(i-1)));
            Vn_new(i) = Pi_d*exp(z(i)) - f + beta*(1-delta)*...
                (q*(px*Vx(i+1)+(1-px)*Vn(i+1))+(1-q)*Vn(i-1));
        end

        Vx_new(1) = (1+tau^(1-rho))*Pi_d*exp(z(1)) - f - fx - fxs +...
            beta*(1-delta)*(q*(Vx(2) + fxs) + (1-q)*0);
        Vn_new(1) = Pi_d*exp(z(1)) - f + beta*(1-delta)*...
                (q*(px*Vx(2)+(1-px)*Vn(2))+(1-q)*0);
        
        Vx_new(N) = (1+tau^(1-rho))*Pi_d*exp(z(N)) - f - fx - fxs +...
            beta*(1-delta)*(q*(Vx(N) + fxs) + (1-q)*((1-pn)*(Vx(N-1)+...
            fxs)+pn*Vn(N-1)));
        Vn_new(N) = Pi_d*exp(z(N)) - f + beta*(1-delta)*...
                (q*(px*Vx(N)+(1-px)*Vn(N))+(1-q)*Vn(N-1));

        Vn_new(Vn_new < 0) = 0;
        Vx_new(Vx_new < 0) = 0;

        abs1 = max(abs(Vn_new - Vn));
        abs2 = max(abs(Vx_new - Vx));
        Vn = Vn_new;
        Vx = Vx_new;
    end


    for i=1:N
        if z(i)< z_x
            EVN(i) = Vn(i)*G(i+1);
            EVX(i) = 0;
        else
            EVN(i) = 0;
            EVX(i) = Vx(i)*G(i+1); 
        end
    end
    EVX_sum = sum(EVX);
    EVN_sum = sum(EVN);
        
    
    if beta*(EVN_sum + EVX_sum) - fe > 0
        Pi_d_max = Pi_d; 
    else
        Pi_d_min = Pi_d;
    end

    Pi_d = (Pi_d_max + Pi_d_min)/2;
    err = Pi_d_max - Pi_d_min;
    iter = iter+1;

end

%% Mass of firms
exit = zeros(N,1);
exit(Vn>0) = 1;
z_exit = find(exit,1);


err1 = 1;
iter1 = 1;
Mn = linspace(0,0,N);
Mx = linspace(0,0,N);
Mn_new = linspace(0,0,N);
Mx_new = linspace(0,0,N);
while (err1 >= tol) && (iter1 <= max_iter)

    for i=1:N
        ln(i) = (rho-1)*Pi_d*exp(z(i+1))*Mn(i);
        lx(i) = (rho-1)*Pi_d*(1+tau^(1-rho))*exp(z(i+1))*Mx(i);
    end
    Lp_Me = sum(ln) + sum(lx);
    M_e = 1/(Lp_Me + fe +f*sum(Mn) + (f+fx+fxs)*sum(Mx));

    for i=2:N-1
        if (z(i+1) > z_n) && (z(i+1) ~= z_x)
            Mx_new(i) = M_e*G(i+1) + (1-delta)*(q*Mx(i-1)+(1-q)*Mx(i+1));
        elseif z(i+1) == z_x
            Mx_new(i) = M_e*G(i+1) + (1-delta)*(q*(Mx(i-1)+Mn(i-1))...
                +(1-q)*Mx(i+1));
        else
            Mx_new(i) = 0;
        end

        if (z(i+1) > z(z_exit)) && (z(i+1) < z_x) && (z(i+1) ~= z_n)
            Mn_new(i) = M_e*G(i+1) + (1-delta)*(q*Mn(i-1)+(1-q)*Mn(i+1));
        elseif z(i+1) == z_n
            Mn_new(i) = M_e*G(i+1) + (1-delta)*(q*Mn(i-1)+(1-q)...
                *(Mx(i+1)+Mn(i+1)));
        else
            Mn_new(i) = 0;
        end
    end

    if (z(2) > z_n)
            Mx_new(1) = M_e*G(2) + (1-delta)*(1-q)*Mx(i+1);
    else
            Mx_new(1) = 0;
    end

    if (z(2) > z(z_exit)) && (z(2) < z_x) && (z(2) ~= z_n)
            Mn_new(1) = M_e*G(2) + (1-delta)*(1-q)*Mn(2);
    elseif z(2) == z_n
            Mn_new(1) = M_e*G(2) + (1-delta)*(1-q)*(Mx(2)+Mn(2));
    else
            Mn_new(1) = 0;
    end


    if (z(N+1) > z_n) && (z(N+1) ~= z_x)
        Mx_new(N) = M_e*G(N+1) + (1-delta)*q*Mx(N-1);
    elseif z(N+1) == z_x
        Mx_new(N) = M_e*G(N+1) + (1-delta)*q*(Mx(N-1)+Mn(N-1));
    else
        Mx_new(N) = 0;
    end

    if (z(N+1) > z(z_exit)) && (z(N+1) < z_x)
        Mn_new(N) = M_e*G(N+1) + (1-delta)*q*Mn(N-1);
    else
        Mn_new(N) = 0;
    end


        

    err1 = max(max(abs(Mn_new-Mn)), max(abs(Mx_new-Mx)));
    iter1 = iter1+1;
    Mn = Mn_new;
    Mx = Mx_new;
end


plot(z(2:N+1),log(ln));
hold on;

plot(z(2:N+1),log(lx));
hold off;
