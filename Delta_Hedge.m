%We now tend to the Delta hedge for the short position in the
%TwinWin Certificate. Assume an amount of N is invested by a client
N=1000000;

%We place the discounted amount on a risk free bankaccount
bank = exp(-r*T)*N;

%In the case that the stock price drops below H, the payout to the investor
%is equal to the stock price at maturity.

%Delta approximation
h = [0.001, 0.01, 0.02, 0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1];
dS_n = size(h,2);
Delta =zeros(1,dS_n);

K_DAO = ceil(S0*exp(r*T));%=65
H_DAO = 60;
DAO_Put_Payoff = max(0,K_DAO-S(:,end)).*max(0,(min(S,[],2)-H_DAO)./abs(min(S,[],2)-H_DAO));

for i=1:dS_n
    S_0 = S0 + h(:,i);
    S1 = zeros(m, n+1);
    v1 = zeros(m, n+1);
    S1(:,1) = S_0;
    v1(:,1) = sigma0^2;
    for j=2:n+1
        S1(:,j)=S1(:,j-1).*(1+(r-q).*dt+sqrt(v1(:,j-1).*dt).*eps1(:,j-1));
        v1(:,j)=abs(v1(:,j-1)+(kappa.*(eta-v1(:,j-1))-theta.^2/4).*dt+theta.*sqrt(v1(:,j-1).*dt).*eps2(:,j-1)+(theta.^2/4).*dt.*eps2(:,j-1).^2);
    end 
    Delta(:,i) = (exp(-r*T)*mean(max(0,S1(:,end))+2*max(0,K_DAO-S(:,end)).*max(0,(min(S1,[],2)-H_DAO)./abs(min(S1,[],2)-H_DAO)))-exp(-r*T)*mean(max(0,S(:,end))+2.*DAO_Put_Payoff))./h(:,i);
end
close all;
figure(1)
hold on
plot(h,floor(N/S0)*Delta,'.','MarkerSize', 10)
xlabel('h')
ylabel('Delta')
hold off