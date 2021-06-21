function main_eurooption
tic
%----- model parameters for diffusion
% dS/S = rate dt + vol (S/spotS)^beta dW
rate = 0.1; % risk free rate
sigma_0 = 0.25; % diffusion vol
theta = 0.2;
k = 10;
beta = -2;
r_0 = 0.01;
r_1 = 0.09;
a_0 = -1;
spotS = 100; % starting level
%----- state-space parameters
N = 200; % number of points in the state-space of X
minS = 1; % smallest element in the grid
maxS = 190; % largest element in the grid
%----- contract details
%T_0 = 1; % bond maturity
T = 1;  % option maturity
strikes = spotS * 0.95; % vector of strikes
m=10;
%B=5;
L = 90; % lower barrier
U = 120; % upper barrier

%--------------------------------------------------------------
% Markov generator algorithm
%--------------------------------------------------------------
% find the grid for the underlying
s_grid1 	= 	generate_grid_for_double_barrier(N, minS, L, spotS, U, maxS, ...
    1, 100, ...
    10, 10, ...
    100, 1);
% set up the call payoff
for y = 1:size(strikes,2)
    payoff(:,y) = max(s_grid1 - strikes(y),0)';
end
%payoff(:,strikes) = max(s_grid1 - strikes,0)';
idx = find(minS<=s_grid1 & s_grid1<=maxS);
newspotS_index = find(s_grid1 == spotS);
barP1 =0;
for l = 1:m
    Tl_1 = T/m*(l-1); %partition
    MG = generate_diffusion_MG1(s_grid1,spotS,beta,a_0,r_0,r_1,sigma_0,theta,k,rate,Tl_1);
    A = T/m*MG;
    %BarrierMG 	= 	MG(idx,idx);
    %if Tl<T
    %    A = T/m*MG;
    %elseif  (Tl-T)>=0 & (Tl-T)<=T_0/m
    %    A = (T_0/m-(Tl-T))*MG;
    %end
   barP1 =barP1+A;
end
A=expm(barP1);
barP = A (newspotS_index,:);
%barPayoff = payoff(idx,:);
barPayoff = payoff(idx,:);
euroPrices = exp(-rate * T) * barP * barPayoff;
 toc 	 
fprintf('\n') 	 
fprintf('---------------- Call 	Options 	-----------------\n') 	 
fprintf('time 	to 	maturity:') 	 
disp(T); 	 
fprintf('strikes:') 	 
disp(strikes); 	 	 
fprintf('call option prices:') 	 
disp(euroPrices);