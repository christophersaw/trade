cd /Users/christophersaw/Desktop/tradeps3;
clc;
rng('default');
T=100;                          % No. of time periods
NJ=200;                         % No. of US states (50) X sectors (4)
u_dot=ones(NJ,T+2);             % Initial guess of u_dot
u_dot_update=zeros(NJ,T+2);
mu=zeros(NJ,NJ,T+1);
mu(:,:,1)=mu1;
L=zeros(NJ,T+1);
L(:,1)=reshape(L0,[NJ,1]);
L_dot=ones(NJ,T)
xi=xi(1,[1,87]);

global X0 pi0 L0 mu1 gamma_nj gamma_njnk_res iota xi alpha

% PREPARE MATRICES

% pi0 is 348 x 87 (NJ by N)
% rows: demand, columns: supply
% I want to know the flow from country n to country i for a given k
pi0_res=zeros(4,87,87); %(sector,supply, demand)
for i=0:86
    for n=1:87
        for j=1:4
            pi0_res(j,n,i+1)=pi0(4*i+j,n); % note the switch from col supply to row supply
        end
    end
end

% gamma_njnk is 348 by 4 (NJ by J)
% rows: supply, columns: demand
% I want to know the flow from sector j to sector k for a given n
gamma_njnk_res=zeros(4,4,87); %(supply, demand, region)
for j=1:4
    for k=1:4
        for n=0:86
            gamma_njnk_res(j,k,n+1)=gamma_njnk(4*n+j,k);
        end
    end
end

% X0 is 348 by 1, reshape for USA (4 sectors, 50 states)
X0_res=zeros(4,50);
for n=0:49
    for j=1:4
        X0_res(j,n+1)=X0(4*n+j);
    end
end

% FIND INTIAL WAGE (use X0 to find w0)
w0=ones(4,50); % guess
w0=fsolve( @(x)initialwage(x),w0);

% PREPARE ITERATION

tol=0.001;
diff=1;

while diff > tol
    % SEQUENTIAL EQUILIBRIUM (OUTER LOOP)
    % Use initial guess of u_dot and solve for mu_t+1 with equation 16
    for t = 1:T
        for y = 1:NJ
                mu(:,y,t+1)=( mu(:,y,t).*( u_dot(y,t+2)^(beta/nu) ) ) / sum( mu(:,y,t).*( u_dot(y,t+2)^(beta/nu) ) );
        end
    end
    
    % Calculate L_t+1 with equation 18
    for t = 1:T
        for y = 1:NJ
                L(y,t+1)=sum( mu(:,y,t) .* L(:,t) );
        end
    L_dot(:,t) = L(:,t+1)/L(:,t);
    end
    
    % Reshape L_dot from (NJ,T) to (J,N,T)
    L_dot_res=zeros(4,50,T);
    for n=0:49
        for j=1:4
            for t=1:T
                L_dot_res(j,n+1,t)=L_dot(4*n+j,t);
            end
        end
    end
    
    
    % TEMPORARY EQUILIBRIUM (INNER LOOP) - incomplete
    % Initiate guess for w_dot 
    w_dot_res=ones(4,50,T);          % index t is t+1 in CDP
    w_dot_res=fsolve( @(x)initialwage(x),w_dot_res);

    % Build arrays
    P_dot_res=zeros(4,50,T);         % index t is t+1 in CDP
    x_dot_res=zeros(4,50,T);         % index t is t+1 in CDP
    A_dot_res=zeros(4,50,T);         % index t is t+1 in CDP
    pi_seq_res=zeros(4,50,50,T+1);   % T+1 because of initial pi0
    X_seq_res=zeros(4,50,T+1);       % T+1 because of initial X0
    w_seq_res=zeros(4,50,T+1);       % T+1 because of initial w0 (pinned by X0)
    L_seq_res=zeros(4,50,T+1);       % T+1 because of initial L0
    
    % Set intial values in period 0
    pi_seq_res([1:4],[1:50],[1:50],1)=pi0_res([1:4],[1:50],[1:50]);
    X_seq_res([1:4],[1:50],1)=X0_res([1:4],[1:50]);
    w_seq_res([1:4],[1:50],1)=w0([1:4],[1:50]); 
    L_seq_res([1:4],[1:50],1)=L0([1:4],[1:50]);
    
    % Given the guess of w_dot and initial w0, fill out the wage sequence
    for t=1:T 
        w_seq_res(:,:,t+1)=w_dot_res(:,:,t).*w_seq_res(:,:,t);
    end
    
    % Equation 11 
    % the exp() is to revert the log(), cleaner to code a simple sum
    for j=1:4
        for n=1:50
            for t=1:T
                x_dot_res(j,n,t)=exp(  gamma_nj(j,n)*xi(j,n)*log(L_dot_res(j,n,t)) ...
                                     + gamma_nj(j,n)*log(w_dot_res(j,n,t) +        ...
                                                         sum ( gamma_njnk_res(:,j,n).*log(P_dot_res(:,n,t)) )...
                                                         ) ...
                                    );
            end
        end
    end
    
    % Equation 12
    for j=1:4
        for n=1:50
            for t=1:T
                P_dot_res(j,n,t)=(...
                                   sum(  pi0_res(:,n,j).*(x_dot_res(j,:,t).^(-theta(j))).*...
                                         ( A_dot_res(j,:,t).^( theta(j).*gamma(j,:) ) )   ...
                                       )...
                                 )^(-1/theta(j));
            end
        end
    end
    
    % Equation 13
    for n=1:50
        for i=1:50
            for t=1:T
                pi_seq_res(:,n,i,t+1)=pi_seq_res(:,n,i,t).*...
                                      ( x_dot_res(:,i,t) ./ P_dot_res(:,n,t) ).^( -theta(:) ).*...
                                      A_dot_res(:,i,t).^( theta(:).*gamma(:,i) ) ;
            end
        end
    end
    
    % Equation 14 - incomplete
    chi_vec=zeros(50,T);
    chi=zeros(T);          % index t is t+1 in CDP
    for i=1:N
        for t=1:T
            chi_vec(i,t)= xi(i)./(1-xi(i)) .* sum(  w_dot_res(:,i,t).*L_dot_res(:,i,t).*w_seq_res(:,i,t).*L_seq_res(:,i,t)  ); 
        end
    end
    for t=1:T
        chi(t)=sum(chi_vec(:,t));
    end
    
    % Want to find X_t+1 given w_t,L_t, w_dot, L_dot
    % This inversion is non-trivial...
    
    % Equation 15
    LHS=zeros(4,50,T);
    RHS=zeros(4,50,T);
    LHS(:,:,:)=w_dot_res(:,:,:).*L_dot_res(:,:,:).*w_seq_res(:,:,:).*L_seq_res(:,:,:);
    for j=1:4
        for n=1:50
            for t=1:T
                RHS(j,n,t)=gamma_nj(j,n)*(1-xi(n))*sum_piX(j,n,t+1);
            end
        end
    end
    Z0=LHS-RHS;

    % Update u_dot with equation 17
    for t = 1:T
        for y = 1:NJ
        omega_dot(:,t+1)=w_dot(:,t+1)/P_dot(:,t+1)
        u_dot_update(:,t+1)=omega_dot(:,t+1)* ( sum( mu(:,y,t).*( u_dot(y,t+2,h)^(beta/nu) ) ) )^nu
        end
    end
    diff=norm(u_dot_update-u_dot)
    u_dot=u_dot_update
end