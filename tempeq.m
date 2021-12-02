function [Z0]=tempeq(w_dot_res,L_dot_res,w0,L0,X0_res,pi0_res)
global gamma_nj gamma_njnk_res iota xi alpha
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

% Equation 14
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

% sum_piX=zeros(4,50,T+1);
% for n=1:50
%     for k=1:4
%         for t=1:T+1
%             sum_piX(k,n,t)=sum(pi_seq_res(k,n,:,t).*X_seq_res(k,:,t));
%         end
%     end
% end
% 
% for n=1:50
%     for j=1:4
%         for t=1:T
%             X_seq_res(j,n,t+1)= sum( gamma_njnk_res(j,:,n) .* sum_piX(n,:,t+1) )...
%                                  + alpha(j,n)*( sum( w_dot_res(:,n,t).*L_dot_res(:,n,t).*w_seq_res(:,n,t).*L_seq_res(:,n,t) )...
%                                  + iota(n).*chi(t) );
%         end
%     end
% end

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

