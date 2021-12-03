%%% QUESTION 3 AND 4

clc; close all; clear;

tic

ind_counterfactual = 1;
% =1 Baseline economy with constant fundamentals
% =2 All fundamentals are constant and equal to their 2000 values except rising Chinese productivity

ind_initial = 1;
% =1 Use the initial value in CDP paper.
% =2 Use zeros or ones matrix as initial value

J=4;   %number of sectors
N=50;  %number of regions in the world
R=50;  %number of U.S. regions
T=200; %number of quarters

load('base_year_4sectors.mat')

Hvect=ones(N*(J),T); %Initial guess

%%%Dynamic problem%%%
tol_Y=0.001;   %tolerance
maxiter=1E+20; %maximum number of iterations
Y_max=1; 
iter=1;

while (iter <= maxiter) && (Y_max > tol_Y)
    J=4;
    N=50;
    %Solving for the path of migration flows mu
    for t = 1:T
        H(:,:,t) = reshape(Hvect(:,t),J,N);
    end
    
    for t=1:T
        Haux0(:,:,t)=reshape(H(:,:,t),(J)*N,1);
    end
    Haux= permute(Haux0,[2,1,3]);
    
    Haux2=zeros((J)*N,N*(J),T);
    for t=1:T
        Haux2(:,:,t)=repmat(Haux(:,:,t),N*(J),1);
    end
    
    %Computing mu1
    num=mu1.*(Haux2(:,:,2).^beta);
    den=sum(num')';
    den=den*ones(1,N*(J));
    mu00=num./den;
    mu=zeros(N*(J),N*(J),T);
    mu(:,:,1)=mu00;  
    
    %Solving for mu(t+1)
    for t=1:T-2
        num=mu(:,:,t).*(Haux2(:,:,t+2).^beta);
        den=sum(num')';
        den=den*ones(1,N*(J));
        mu(:,:,t+1)=num./den;  
    end
    
    %Solving for the path of employment
    L00=reshape(L0,J,N);
    L00_aux2=reshape(L00,N*(J),1);
    L00_aux4=repmat(L00_aux2,1,N*(J));
    L00_aux5=mu00.*L00_aux4;
    L1=sum(L00_aux5)';
    L1=reshape(L1,J,N);
    Ldyn=zeros(J,N,T);
    Ldyn(:,:,2)=L1;
    Ldyn(:,:,1)=reshape(L0,(J),N);
    for t=2:T-1
        aux=reshape(Ldyn(:,:,t),N*(J),1);
        aux3=repmat(aux,1,N*(J));
        aux4=mu(:,:,t).*aux3;
        aux5=sum(aux4)';
        Ldyn(:,:,t+1)=reshape(aux5,J,N);   
    end
    Ldyn(:,:,T)=0;
    
    %%%%Temporary Equilibrium%%%%%
    
    J=4;
    N=87;
    vfactor=-0.05;  %adjustment factor of the wage update
    tol=0.001;      %tolerance
    maxiter=1E+20;  %maximum number of iterations
    
    if ind_initial == 1
        load('initial1.mat')
    else
        load('initial2.mat')
    end
    
    realwages=ones(J,N,T);  %real wages. This matrix will store equilibrium real wages from the temporary equilibrium at each time t
    
    A_dot=ones(J,N,T); % fundamental TFP
    if ind_counterfactual == 2
        n_china = 57;
        period = 31;
        A_end = 44;
        
        A_dot(1,n_china,1:period) = A_end^(1/period);
    end
    
    %static sub-problem at each time t
    for t=1:T-2
        disp(t);
        %Shocks
        kappa_hat=ones(J*N,N);                      % relative chane in trade costs
        A_hat=A_dot(:,:,t);                         % relative change in technology
        Snp = zeros(N,1);                           % trade balance at t=t+1
        
        % Initialize vectors of factor prices (om) and good prices (p)
        om=ones(J,N);                               % initial guess for factor prices
        pf0=ones(J,N);                              % initial guess for good prices
        Ljn_hat=ones(J,N);                          % change in employment
        Ljn_hat(:,1:R)=Ldyn(:,:,t+1)./Ldyn(:,:,t);  % change of employment in the US
        
        ommax = 1; itw = 1;
        while (itw <= maxiter) && (ommax > tol)
            
            %%%Calculating good prices and input bundle consistent with factor prices%%%
            pf0=ones(J,N);
            
            pfmax=1;
            
            it       = 1;
            while (it <= maxiter) && (pfmax > tol)
                lom=log(om);
                lp=log(pf0);
                
                %Calculating input bundle costs
                for i=1:N
                    lc(:,i)=gamma_nj(:,i).*lom(:,i)+(gamma_njnk(1+(i-1)*J:J*i,:)'*lp(:,i));
                end
                c=exp(lc);
                for   j=1:1:J
                    idx=1+(j-1)*N:1:N*j;
                    LT(idx,1)=ones(N,1)*theta(j);
                end
                Din_k=pi0.*(kappa_hat.^(-1./(LT*ones(1,N))));
                
                %Calculating change in prices
                for j=1:1:J
                    for n=1:1:N
                        phat(j,n)=Din_k(n+(j-1)*N,:)*((A_hat(j,:).^(gamma_nj(j,:)./(theta(j)))).*(c(j,:).^(-1/theta(j))))';
                        phat(j,n)=phat(j,n)^(-theta(j));
                    end
                end
                
                pfdev    = abs(phat - pf0); %checking tolerance
                pf0      = phat;
                pfmax    = max(max(pfdev));
                it       = it + 1;
            end
            
            %%%Calculating bilateral trade shares%%%
            
            %Reformating theta vector
            for   j=1:1:J
                idx=1+(j-1)*N:1:N*j;
                LT(idx,1)=ones(N,1)*theta(j);
            end
            
            %Calculating bilateral trade shares
            for n=1:1:N
                cp(:,n)=c(:,n).^( - 1./theta );
                phatp(:,n)=phat(:,n).^(-1./theta);
            end
            
            Din_k = pi0.*(kappa_hat.^(-1./(LT*ones(1,N))));
            
            for n=1:1:N
                idx=n:N:length(pi0)-(N-n);
                DD(idx,:) = Din_k(idx,:).*(cp.*(A_hat.^(gamma_nj./(theta*ones(1,N)))));
            end
            
            for n=1:1:N
                idx=n:N:length(pi0)-(N-n);
                pi(idx,:)=DD(idx,:)./(phatp(:,n)*ones(1,N));
            end
            
            %%%Calculating total expenditure%%%
            
            VARjnp = (VARjn0.*om .*(Ljn_hat.^(1-xi)));  % Structures in the counterfactual equilibrium
            VARp=sum(VARjnp)';
            Chip = sum(VARp) ;                          % Total revenues in the global portfolio
            Bnp = Snp - iota.*Chip + VARp;              % New trade balance
            
            % Calculating the Omega matrix
            
            NBP = zeros(size(pi'));
            
            for j = 1:1:N
                for n = 1:1:N
                    NBP(j , 1 + (n-1)*J : n*J) = pi([n:N:J*N], j);
                end
            end
            
            NNBP=kron(NBP,ones(J,1));
            GG=kron(ones(1,N),gamma_njnk);
            GP=GG.*NNBP;
            OM=eye(J*N,J*N)-GP ;
            
            % Calculating total expenditures
            
            aux=( sum(om.*(Ljn_hat.^(1-xi)).*(VARjn0 + VALjn0))' -Bnp);
            aux2=kron(aux,ones(J,1));
            X=inv(OM)*(reshape(alpha,N*J,1).*aux2);
            Xp=reshape(X,J,N);
            
            %%%Calculating new wages using the factor market clearing condition%%%
            
            %Calculating new factor prices using the labor market clearing condition
            PQ_vec=reshape(Xp',1,J*N)';
            for n=1:1:N
                DP(:,n)=pi(:,n).*PQ_vec;
            end
            for j=1:1:J
                for n=1:1:N
                    Exjnp(j,n)=sum(DP(1+N*(j-1):N*j,n))';
                end
            end
            
            aux4=gamma_nj.*Exjnp;
            aux5=(aux4);
            omef0=ones(J,N);
            omef0(:,1:R) = aux5(:,1:R)./((Ljn_hat(:,1:R).^(1-xi(:,1:R))).*(VARjn0(:,1:R) + VALjn0(:,1:R)));
            VAR=sum(VARjn0)';
            VAL=sum(VALjn0)';
            aux5=sum(aux4)';
            omef0(:,R+1:N)=ones(J,1)*(aux5(R+1:N)./(VAR(R+1:N) + VAL(R+1:N)))';
            % Excess function
            ZW = (om-omef0);
            % Iteration factor prices
            om1 =om.*(1+vfactor*ZW./om);
            om11=reshape(om1,N*J,1);
            om00=reshape(om,N*J,1);
            omusa=om11(1:J*R,1)-om00(1:J*R,1);
            omrow=om1(1,R+1:N)'-om(1,R+1:N)';
            omworld=[omusa; omrow];
            ommax=sum(abs(omworld).^2); %checking tolerance
            om=om1;
            itw=itw+1;
        end
        
        % Recovering wages and rental rates
        wf0=zeros(J,N);
        wf0(:,1:R)=om(:,1:R).*(Ljn_hat(:,1:R).^(-xi(:,1:R)));
        wf0(:,R+1:N)=om(:,R+1:N);
        
        % Recovering deficits
        VARjnp = (VARjn0.*om .*(Ljn_hat.^(1-xi)));
        
        % Recovering other variables
        VALjn0 = (wf0.*Ljn_hat.*VALjn0);
        
        Phat=prod(phat.^(alpha)); %price index
        
        VARjn0 = VARjnp;
        pi0 = pi;
        
        %%%%updating the initial conditions%%%
        
        %storing equilibrium real wages
        realwages(:,:,t+1)=wf0./(ones(J,1)*Phat);
        
    end
    
    %%%%Solving for the new path of values in time differences%%%
    realwagesaux=zeros(J,N,T);
    for t=1:T
        realwagesaux(1,:,t)=1;
        realwagesaux(:,:,t)=realwages(:,:,t);
    end
    
    rwage=ones(R,J,T);
    
    for t=1:T
        rwage(:,:,t)=permute(realwagesaux(:,1:R,t),[2,1,3]);
    end
    
    rw=rwage.^(1/nu);
    rw=permute(rw,[2,1,3]);
    rw_aux=zeros(R*(J),1,T);
    for t=1:T
        rw_aux(:,:,t)=reshape(rw(:,:,t),R*(J),1);
    end
    
    rw_aux2=zeros(R*(J),R*(J),T);
    for t=1:T
        rw_aux2(:,:,t)=repmat(rw_aux(:,:,t),1,R*(J));
    end
    
    rwagenu=zeros(R*(J),R*(J),T);
    
    for t=1:T-1
        rwagenu(:,:,t+1)=mu(:,:,t).*rw_aux2(:,:,t+1);
    end
    
    num=zeros(size(rwagenu));
    for t=1:T-1
        num(:,:,t)=rwagenu(:,:,t).*(Haux2(:,:,t+1).^beta);
    end
    
    Y=zeros(R*(J),1,T);
    for t=1:T
        Y(:,:,t)=sum(num(:,:,t)')';  
    end
    Y(:,:,T)=1;
    
    Ynew=zeros((J)*R, T);
    for t=1:T
        Ynew(:,t)=Y(:,:,t);
    end
    Ynew(:,T)=1;
    Hvect(:,T)=1;
    
    %Excess function
    check=zeros (T,1);
    for t=2:T
        checkY(t,1)=max(abs(Ynew(:,t)-Hvect(:,t)));
    end
    Y_max=max(checkY)
    Ymax0=Y_max;
    Hvect=0.5*Ynew+0.5*Hvect;
    iter=iter+1;
end
toc

if ind_counterfactual == 1
    save('result1.mat', 'Ynew', 'Ldyn');
else
    save('result2.mat', 'Ynew', 'Ldyn');
end

%%% QUESTION 5

clc; close all; clear;

T = 200;

load("result1.mat")

Ldyn1 = Ldyn;

load("result2.mat")

Ldyn2 = Ldyn;

manu_emp_change = (sum(Ldyn2(1,:,T-1)) - sum(Ldyn1(1,:,T-1))) / sum(Ldyn1(1,:,T-1));