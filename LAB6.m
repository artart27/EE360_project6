%% Lab 5
%
%   Name: Paco Ellaga            (SID: 009602331)
%   Name: Arturo Ramos Hernandez (SID:017264297)
%   Date: November 2, 2018
%
clear, clc
%% Part 1
% Create the plot of Figure 1, showing the effect on sample size on the 
% confidence intervals.Use the following parameters:
clear
close all
Nb = 1000000;        % Total number of bearings:
mu_p = 75;           % population mean in grams
sig_p  = 7.5;        % population standard deviation in grams, sigma
n = 1:200; n = n';   % sample sizes: n = 1,2,...,200

kmax = max(size(n));
M = zeros(length(n),2);
Bearings = sig_p*randn(Nb,1)+mu_p;
for k = 1:kmax
        x_index=ceil(Nb*rand(k,1));
        x = Bearings(x_index);
        x_bar=mean(x);
        M(k,:)=[k x_bar ];
  
end

% plot showing the effect on samplesize on the confidence intervals.
figure(1); %plot(n,M(:,2),'o', n, mu_p*ones(size(n)),'k')
figure(1); hold on
% 1.96 corresponds to the 95% confidence interval
% 2.58 corresponds to the 95% confidence interval
% upper bound  of confidence interval
plot(n, mu_p+1.96*sig_p./sqrt(n),'r', n, mu_p+2.58*sig_p./sqrt(n),'m--')
% lower bound of confidence interval
plot(n, mu_p-1.96*sig_p./sqrt(n),'r', n, mu_p-2.58*sig_p./sqrt(n),'m--')
figure(1); hold off

%% Part 2
%
% n = 5
%
n = 5;
A1 = 0;
B1 = 0;
C1 = 0;
D1 = 0;
Z = 0;
while Z < 10000

k = randi(1e6 - n);
sample = Bearings(k:k+n);
mu = (sum(sample))/n;
sigma = sum(((sample-mu).^2))/(n-1);
% Normal distribuiton at 95% confidence
mu_n_low_95 = mu-1.96*(sigma/sqrt(n));
mu_n_up_95 = mu+1.96*(sigma/sqrt(n));

if (mu_n_low_95<mu_p)&&(mu_n_up_95>mu_p)
    A1 = A1 + 1;
end

% t distribuiton at 95% confidence
mu_t_low_95 = mu-2.78*(sigma/sqrt(n));
mu_t_up_95 = mu+2.78*(sigma/sqrt(n));

if (mu_t_low_95<mu_p)&&(mu_t_up_95>mu_p)
    B1 = B1 + 1;
end

% Normal distribuiton at 99% confidence
mu_n_low_99 = mu-2.58*(sigma/sqrt(n));
mu_n_up_99 = mu+2.58*(sigma/sqrt(n));

if (mu_n_low_99<mu_p)&&(mu_n_up_99>mu_p)
    C1 = C1 + 1;
end
% t distribuiton at 99% confidence
mu_t_low_99 = mu-4.60*(sigma/sqrt(n));
mu_t_up_99 = mu+4.60*(sigma/sqrt(n));

if (mu_t_low_99<mu_p)&&(mu_t_up_99>mu_p)
    D1 = D1 + 1;
end

Z = Z + 1;
end

%%
% n = 40; t at 95: 2.02; t at 99: 2.70;
%
n = 40;
k = randi(1e6 - n);
sample = Bearings(k:k+n);
mu = (sum(sample))/n;
sigma = sum(((sample-mu).^2))/(n-1);
% Normal distribuiton at 95% confidence
mu_normal_low_95 = mu-1.96*(sigma/sqrt(n));
mu_normal_up_95 = mu+1.96*(sigma/sqrt(n));
% t distribuiton at 95% confidence
mu_t_low_95 = mu-2.02*(sigma/sqrt(n));
mu_t_up_95 = mu+2.02*(sigma/sqrt(n));

% Normal distribuiton at 99% confidence
mu_normal_low_99_n5 = mu-2.58*(sigma/sqrt(n));
mu_normal_up_99_n5 = mu+2.58*(sigma/sqrt(n));
% t distribuiton at 99% confidence
mu_t_low_99_n5 = mu-2.70*(sigma/sqrt(n));
mu_t_up_99_n5 = mu+2.70*(sigma/sqrt(n));

%%
% n = 120; t at 95: 1.98; t at 99: 2.62;
%
n = 120;
k = randi(1e6 - n);
sample = Bearings(k:k+n);
mu = (sum(sample))/n;
sigma = sum(((sample-mu).^2))/(n-1);
% Normal distribuiton at 95% confidence
mu_normal_low_95 = mu-1.96*(sigma/sqrt(n));
mu_normal_up_95 = mu+1.96*(sigma/sqrt(n));
% t distribuiton at 95% confidence
mu_t_low_95 = mu-1.98*(sigma/sqrt(n));
mu_t_up_95 = mu+1.98*(sigma/sqrt(n));

% Normal distribuiton at 99% confidence
mu_normal_low_99_n5 = mu-2.58*(sigma/sqrt(n));
mu_normal_up_99_n5 = mu+2.58*(sigma/sqrt(n));
% t distribuiton at 99% confidence
mu_t_low_99_n5 = mu-2.62*(sigma/sqrt(n));
mu_t_up_99_n5 = mu+2.62*(sigma/sqrt(n));
%%
% Result tabulate
% 
% Sample_size = [5;40;120];
% N_dis_95= [n_5_95;n_40_95;n_120_95];
% N_dis_99= [n_5_99;n_40_99;n_120_99];
% t_dis_95= [t_5_95;t_40_95;t_120_95];
% t_dis_99= [t_5_99;t_40_99;t_120_99];
% T = table(Sample_size,N_dis_95,N_dis_99,t_dis_95,t_dis_99);
% disp(T)