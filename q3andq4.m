
S0 = 62.96;

n1 = 1000;
n2 = 2000;
n3 = 1000;

mu = 0.0516;
sig = 0.1988;
r = 0.03;
K = 64;
T = 14/52;
t = 1/52;
N_Monte = 10000;
N_Hist = 261;

%3a 
%ii
delta_A = delta_call(S0, K, r, sig, 0, T)
theta_A = theta_call(S0, K, r, sig, 0, T)
delta_B = delta_put(S0, K, r, sig, 0, T)
theta_B = theta_put(S0, K, r, sig, 0, T)

%iv
VaR95 = VaRalpha(0.95, mu, sig, S0, t)
VaR99 = VaRalpha(0.99, mu, sig, S0, t)
ES95 = ESalpha(0.95, mu, sig, S0, t)
ES99 = ESalpha(0.99, mu, sig, S0, t)

alpha = 0:0.01:1;
figure;
plot(alpha, VaRalpha(alpha, mu, sig, S0, t), '-r', 'Linewidth', 2);
title('VaR vs \alpha');
xlabel('\alpha');
ylabel('VaR');

%v
gamma_AB = gamma_option(S0, K, r, sig, 0, T)




%3b

D = readmatrix('BNS.TO.weekly.xlsx','Range','B2:B263');
Retrn=zeros(261,1);
for i=1:261
    Retrn(i)=D(i+1)/D(i);
end
HistS = S0.*Retrn;
HistVt = n1.*HistS + n2.*callOptA(HistS, t, sig, r, K, T) + n3.*putOptB(HistS,t,sig,r,K,T);
HistV0 = Vt(n1, n2, n3, mu, sig, 0, S0, r, K, T, N_Hist);
L_Hist = -(HistVt-HistV0);

%--i
figure;
histogram(L_Hist, 'Normalization', 'probability');
title('Historical Simulated Losses');
xlabel('L');
ylabel('Relative Frequency');

%--ii
VaR95_Hist = quantile(L_Hist, 0.95) 
ES95_Hist = mean(L_Hist(L_Hist>VaR95_Hist))

%--iii
L_Linearized_Hist = linearizedLossHist(n1, n2, n3, sig, t, HistS, S0, r, K, T);

figure;
histogram(L_Linearized_Hist, 'Normalization', 'probability');
title('Historical Simulated Linearized Losses');
xlabel('L Linearized');
ylabel('Relative Frequency');
VaR95_Hist_Linearized = quantile(L_Linearized_Hist, 0.95) 
ES95_Hist_Linearized = mean(L_Linearized_Hist(L_Linearized_Hist>VaR95_Hist_Linearized))

%3c
%i
S = logNormS(S0, mu,sig,t,N_Monte);
vt = Vt(n1, n2, n3, mu, sig, t, S, r, K, T, N_Monte);
v0 = Vt(n1, n2, n3, mu, sig, 0, S0, r, K, T, N_Monte);
L_Monte = lossRND(vt, v0);

figure;
histogram(L_Monte, 'Normalization', 'probability');
title('1-week ahead Loss Histogram');
xlabel('L');
ylabel('Relative Frequency');

VaR95MonteSim = quantile(L_Monte, 0.95)
ES95MonteSim = mean(L_Monte(L_Monte > VaR95MonteSim))

VaR99MonteSim = quantile(L_Monte, 0.99)
ES99MonteSim = mean(L_Monte(L_Monte > VaR99MonteSim))

%--ii

L_deltaGamma = gammaDeltaApproxL(n1, n2, n3, mu, sig, t, S0, r, K, T, N_Monte);
VaR95deltaGamma = quantile(L_deltaGamma, 0.95)
ES95deltaGamma = mean(L_deltaGamma(L_deltaGamma > VaR95deltaGamma))

VaR99deltaGamma = quantile(L_deltaGamma, 0.99)
ES99deltaGamma = mean(L_deltaGamma(L_deltaGamma > VaR99deltaGamma))

% L_Lin = linearizedLoss(n1, n2, n3, mu, sig, t, S0, r, K, T, N_Monte);
% VaR95Lin = quantile(L_Lin, 0.95)
% ES95Lin = mean(L_Lin(L_Lin > VaR95Lin))
% 
% VaR99Lin = quantile(L_Lin, 0.99)
% ES95Lin = mean(L_Lin(L_Lin > VaR99Lin))


%Q4 

%---a
%-----i
mean_hist_L = mean(L_Hist);
var_hist_L = var(L_Hist);
fx_value_hist = normpdf(norminv(0.95,mean_hist_L,sqrt(var_hist_L)), mean_hist_L, sqrt(var_hist_L));
st_error_hist = (1/fx_value_hist)*sqrt(0.95*0.005/N_Hist)
lower_CI_VaR95_Hist = VaR95_Hist-1.96*st_error_hist
upper_CI_VaR95_Hist = VaR95_Hist+1.96*st_error_hist

%-----ii
mean_Monte_L = mean(L_Monte);
var_Monte_L = var(L_Monte);
fx_value_Monte = normpdf(norminv(0.95,mean_Monte_L,sqrt(var_Monte_L)), mean_Monte_L, sqrt(var_Monte_L));
st_error_Monte = (1/fx_value_Monte)*sqrt(0.95*0.005/N_Monte)
lower_CI_VaR95_Monte = VaR95MonteSim-1.96*st_error_Monte
upper_CI_VaR95_Monte = VaR95MonteSim+1.96*st_error_Monte

%----4b

S = zeros(N_Monte, 30);
vt = zeros(N_Monte, 30);
v0 = zeros(N_Monte, 30);
L_Monte = zeros(N_Monte, 30);

for i = 1:1:30
    S(:,i) = logNormS(S0, mu,sig,t,N_Monte);
    vt(:,i) = Vt(n1, n2, n3, mu, sig, t, S(:,i), r, K, T, N_Monte);
    v0(:,i) = Vt(n1, n2, n3, mu, sig, 0, S0, r, K, T, N_Monte);
    L_Monte(:,i) = lossRND(vt(:,i), v0(:,i));
end

VaR95MonteSim = zeros(1,30);
ES95MonteSim = zeros(1,30);

for i = 1:1:30
    VaR95MonteSim(i) = quantile(L_Monte(:,i), 0.95);
    temp = L_Monte(:,i);
    ES95MonteSim(i) = mean(temp(L_Monte(:,i) > VaR95MonteSim(i)));
end

avg30_VaR95_Monte = mean(VaR95MonteSim)
std30_VaR95_Monte = std(VaR95MonteSim)
avg30_ES95_Monte = mean(ES95MonteSim)
std30_ES95_Monte = std(ES95MonteSim)

lower_CI_VaR95_Monte30 = avg30_VaR95_Monte-1.96*std30_VaR95_Monte/sqrt(30)
upper_CI_VaR95_Monte30 = avg30_VaR95_Monte+1.96*std30_VaR95_Monte/sqrt(30)
lower_CI_ES95_Monte30 = avg30_ES95_Monte-1.96*std30_ES95_Monte/sqrt(30)
upper_CI_ES95_Monte30 = avg30_ES95_Monte+1.96*std30_ES95_Monte/sqrt(30)

function y = VaRalpha(alpha, mu, sig, S0, t)
    temp = (mu-0.5.*sig^2).*t - sig.*sqrt(t).*norminv(alpha);
	y = 92259+15081.*t-1465.4*S0.*exp(temp);
end

function y = ESalpha(alpha, mu, sig, S0, t)
    temp = log((92549.02-VaRalpha(alpha, mu, sig, S0, t))/(1465.4*S0))
    temp1 = normcdf(-sig*sqrt(t)+(1/sig)*(temp-(mu-0.5*sig^2)*t))
    y = (-1465.4*S0)*exp((mu-0.5*sig^2)*t+0.5*t*sig^2)*temp1/(1-alpha) + 99549.02;
end

function y = logNormS(S0, mu,sig,t,N)
	Xt = (mu - sig.^2).*t + sig.*sqrt(t).*normrnd(0,1,[N,1]);
	y = S0.*exp(Xt);
end

function y = d1(S, K, r, sig, t, T)
	y = (log(S./K)+(r+0.5.*sig^2).*(T-t))/(sig.*sqrt(T-t));
end

function y = d2(d1_call, sig, t, T)
	y = d1_call-sig.*sqrt(T-t);
end

function y = callOptA(S, t, sig, r, K, T)
	d1_A = d1(S, K, r, sig, t, T);
	y = S.*normcdf(d1_A)-K.*exp(-r.*(T-t)).*normcdf(d2(d1_A, sig, t, T));
end

function y = putOptB(S,t,sig,r,K,T)
	d1_B = d1(S, K, r, sig, t, T);
	y = K.*exp(-r.*(T-t)).*normcdf(-1.*d2(d1_B, sig, t, T))-S.*normcdf(-1.*d1_B);
end

function v = Vt(n1, n2, n3, mu, sig, t, S, r, K, T, N)
	v = n1.*logNormS(S, mu, sig, t, N) + n2.*callOptA(S,t,sig,r,K,T) + n3.*putOptB(S,t,sig,r,K,T);
end

function L = lossRND(vt, v0)
	L = -1.*(vt-v0);
end

function y = delta_call(S, K, r, sig, t, T)
	y = normcdf(d1(S, K, r, sig, t, T));
end

function y = theta_call(S, K, r, sig, t, T)
	d1_call = d1(S, K, r, sig, t, T);
	y = -1.*S.*sig.*normpdf(d1_call)/(2.*sqrt(T-t)) - r.*K.*exp(-r.*(T-t)).*normcdf(d2(d1_call, sig, t, T));
end

function y = delta_put(S, K, r, sig, t, T)
	y = -1.*normcdf(-1.*d1(S, K, r, sig, t, T));
end

function y = gamma_option(S, K, r, sig, t, T)
    d1_option = d1(S, K, r, sig, t, T);
    y = normpdf(d1_option)/(S*sig*sqrt(T));
end

function y = theta_put(S, K, r, sig, t, T)
	d1_put = d1(S, K, r, sig, t, T);
	y = -1.*S.*sig.*normpdf(d1_put)./(2.*sqrt(T-t)) + r.*K.*exp(-r.*(T-t)).*normcdf(-1.*d2(d1_put, sig, t, T));
end

function y = linearizedLoss(n1, n2, n3, mu, sig, t, S0, r, K, T, N)
    delta_A = delta_call(S0, K, r, sig, 0, T);
	delta_B = delta_put(S0, K, r, sig, 0, T);
	theta_A = theta_call(S0, K, r, sig, 0, T);
	theta_B = theta_put(S0, K, r, sig, 0, T);
	St = logNormS(S0, mu, sig, t, N);
   
    y = (-n1-n2.*delta_A-n3.*delta_B).*(St-S0) + (-n2.*theta_A-n3.*theta_B).*t;
end

function y = linearizedLossHist(n1, n2, n3, sig, t, St, S0, r, K, T)
    delta_A = delta_call(S0, K, r, sig, 0, T);
	delta_B = delta_put(S0, K, r, sig, 0, T);
	theta_A = theta_call(S0, K, r, sig, 0, T);
	theta_B = theta_put(S0, K, r, sig, 0, T);
	
    y = (-n1-n2.*delta_A-n3.*delta_B).*(St-S0) + (-n2.*theta_A-n3.*theta_B).*t;    
    
end

function y = gammaDeltaApproxL(n1, n2, n3, mu, sig, t, S0, r, K, T, N)
    delta_A = delta_call(S0, K, r, sig, 0, T);
	delta_B = delta_put(S0, K, r, sig, 0, T);
    gamma_AB = gamma_option(S0, K, r, sig, 0, T);
    
    St = logNormS(S0, mu, sig, t, N);
    
    y = (-n1-n2.*delta_A-n3.*delta_B).*(St-S0) + (-1500*gamma_AB).*(St-S0).^2;
    %y = 1465.5.*(S0 - St) - 92.1.*(St-S0).^2;
end

    

