function [Z0]=initialwage(w0)
% SOLVE FOR W0 GIVEN XO, pi0, L0, mu1
global X0 pi0 L0 mu1 gamma_nj gamma_njnk_res iota xi alpha
X0_rhs=zeros(4,50);             % use goods market clearing and structures market clearing conditions to fill X0_rhs
piX0=zeros(348,87);
piX0_res=zeros(87,87,4); %(supply, demand, sector)
for i=1:87
    piX0(i,:)=pi0(i,:).*X0(i); 
end
for i=0:86
    for n=1:87
        for j=1:4
            piX0_res(n,i+1,j)=piX0(4*i+j,n); % note the switch from col supply to row supply
        end
    end
end
sum_piX=zeros(87,4);
for n=1:87
    for k=1:4
        sum_piX(n,k)=sum(piX0_res(n,:,k)); % take sum over all buyers
    end
end
chi_country=zeros(87);
for n=1:87
    chi_country(n)=sum( gamma_nj(:,n).*xi(n).*transpose(sum_piX(n,:)) );  % inside sum over sectors (note this term is meaningless)
end
chi=sum(chi_country(:));                                                       % outside sum over countries

for n=1:50
    for j=1:4
        X0_rhs(j,n)= sum ( gamma_njnk_res(j,:,n) .* sum_piX(n,:) ) + alpha(j,n)*( sum( w0(:,n).*L0(:,n) )+ iota(n).*chi );
    end
end
X0_usa=x0([1:200],1);
Z0=X0_usa-X0_rhs;