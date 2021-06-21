function [kappa_l,z] = generate_diffusion_MG1(S,mu,alpha,beta,rate,Tl_1)
%  Generator for the time-inhomo chain that approximates
%  diffusion: dS/S = rate  dt + vol dW

N = size(S, 2); % number of points in the state-space
Ll(1,:) = 0;         % absorbing boundaries
Ll(N,:) = 0;
kl = zeros(N,N);
kl(1,:) = 1; % 1+kappa
kl(N,:) = 1; % 1+kappa
bar_Ll(1,:) = 0;
bar_Ll(N,:) = 0;
for x = 2:N-1
    dSp = S(x+1)-S(x); 
    dSm = S(x-1)-S(x);
    muS = mu*S(x);%((B-log(S(x)))/(T_0-Tl_1)+1/2*vol^2)*S(x)       
    volS = beta^2*(S(x)-alpha*exp(mu*Tl_1))^2; %beta*(S(x)-alpha*exp(mu*Tl_1));  
    A = [dSp dSm; dSp^2 dSm^2];
    b = [muS; volS];
    ret = A\b;
    Ll(x,x+1) = ret(1);
    Ll(x,x-1) = ret(2);
    Ll(x,x) = -(Ll(x,x+1)+Ll(x,x-1));      %generator under P
    syms a c d
    S1 = -Ll(x,x+1)+1/a+d*Ll(x,x+1)*dSp;
    S2 = -Ll(x,x-1)+1/c+d*Ll(x,x-1)*dSm;
    S3 = c*Ll(x,x-1)*dSm+a*Ll(x,x+1)*dSp-rate*S(x);
    [a,c,d]=solve(S1,S2,S3);
    double([a,c,d]);                              
    bar_ret=[a;c]; 
    kl(x,x+1)=bar_ret(1);
    kl(x,x-1)=bar_ret(2);
    kl(x,x) = 1;                                    %1+kappa
    bar_Ll(x,x+1) = kl(x,x+1)* Ll(x,x+1);
    bar_Ll(x,x-1) = kl(x,x-1)* Ll(x,x-1);
    bar_Ll(x,x) = -(bar_Ll(x,x-1)+bar_Ll(x,x+1));   %generator under Q
end

[row,column] = size(kl);
kappa_l = kl-ones(row,column);
z = bar_Ll;                                         