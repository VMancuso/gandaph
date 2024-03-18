function perf = gandalph(nmec, ncloud, kmec, kcloud, ...
    num_users, lambda_user, ...
    rtt_cli_proxy, rtt_proxy_mec, rtt_proxy_cloud, ...
    Preq, Pans, Pcyc, ...
    alpha_step, timeout) 

%gandalph(1, 10, 10, 50, 50, 25, 24.875, 0.079519, 24.543296, 2000, 4000, 5e6, 0.001, 0.075)

syms a b n k t x y s

%nmec: max in service at the MEC
%ncloud: max in service at the Cloud
%kmec: max jobs ammitable at the MEC, including in service
%kcloud: max jobs admittable at the Cloud, including in service
%num_users: number of users;
%lambda_user: traffci offerd by each user (arr/s);

% total offered traffic:
lambda_in = num_users * lambda_user;

%rtt_cli_proxy: RTT between user and backhaul
%rtt_proxy_mec: RTT between backhaul and MEC
%rtt_proxy_cloud: RTT between backhaul and cloud
rtts = [rtt_cli_proxy rtt_proxy_mec rtt_proxy_cloud]; 

%other parameters

test_points = 3;

%Preq: bits per request packet
%Pans: bits per reply packet
%Pcyc: CPU cycles per request


%RACH parameters 
Tmax = 0.5e-3;
Wmax = 0.5e-3;
kmax = 10; %max number of RACH attempts
EBa = 0.010; %backoff time between RACH attempts, in seconds
tau = 0.001; %interval between RACH opportunities
N = 54; %available preambles

%range of probabilities for edge rouiting
alpha_min=0;
%alpha_step: resolution of the search for the optimum (input parameter)
alpha_max=1;
alpha_in = alpha_min : alpha_step : alpha_max;


%This code uses vectors to refer to job sizes and service rate 
%of the differnt considered network elements:
% position 1: RACH access of the user to the base station
% position 2: Data Uplink between user and backhaul
% position 3: Data Downlink between backhaul and user
% position 4: Uplink from backhaul to MEC server
% position 5: Downlink from MEC to backhaul
% position 6: Uplink from backhaul to Cloud server
% position 7: Downlink from Cloud to backhaul
% position 8: MEC server processing
% position 9: Cloud server processing

%vector of job sizes
P = [Preq Preq Pans Preq Pans Preq Pans Pcyc Pcyc];

%vector of capacities
%C(1) is in arrivals/s (request service capacity of the RACH)
%C(2) to C(7) are expressed in bps, becasue they are for network queues
%C(8) and C(9) are CPU cycles/s (computing capacity of servers)
C(1) = N/exp(1)/tau/10; %RACH collision probabilities are negligible at this level 
C(2:7) = [10 25 20 35 100 100] .* (1e6);
C(8:9) = [1 1] .* (1e9);

%vector of service rates (all expressed as arrivals/s)
mu(1) = C(1);
mu(2:9) = C(2:9)./P(2:9);

%vector of latencies (queueing not included)
d = [0 rtts(1)/2 rtts(1)/2 rtts(2)/2 rtts(2)/2 (rtts(3)-rtts(2))/2 (rtts(3)-rtts(2))/2 2 1] /1000;

%baseline delay at MEC and Cloud, no queueing
del_mec = d(1)+d(2)+d(3)+d(4)+d(5)+d(8);
del_cloud = d(1)+d(2)+d(3)+d(4)+d(5)+d(6)+d(7)+d(9);

%application timeout
%timeout: service timeout in seconds 

%total loss function for x arrivals/s and alpha = a
loss = @(x,a) tloss(x,a,mu(8),nmec,kmec,mu(9),ncloud,kcloud);

%Options for latency distributions:
%RACH
r = 1:kmax;
pa = (1-exp(-r)) .* exp(-(1/2)*r.*(r-1));
clear r;
F1 = @(s) (1-exp(-s*(Tmax+Wmax))) ./ ((Tmax+Wmax)*s) .* sum(pa.*(exp(-s*Tmax)./(1+s*EBa)).^([1:kmax]-1));
pdf1 = matlabFunction(ilaplace(F1,s,x));
cdf1 = matlabFunction(ilaplace(F1(s)./s,s,x));

%M/D/1-PS queue
F2 = @(s,a,b) (1-a/b).*(a+s).^2.*exp(-(a+s)/b)./(s.^2+a.*(s+(a+s).*(1-a/b)).*exp(-(a+s)/b));
f2 = matlabFunction(ilaplace(F2(s,a,b),s,x));
pdf2 = @(x,a,b)piecewise(x<0,0,f2(x,a,b));

%M/D/1-FIFO queue
F3 = @(s,a,b) (1-a/b).*s.*exp(-s/b)./(s-a.*(1-exp(-s/b)));
f3 = matlabFunction(ilaplace(F3,s,x));
f3 = @(x,a,b)piecewise(x<0,0,f3(x,a,b));

%M/M/1 queue
F4 = @(s,a,b) (b-a)./(s+(b-a));
pdf4 = @(x,a,b) (b-a) .* exp(-(b-a).*x);
cdf4 = @(x,a,b) (1 - exp(-(b-a).*x));

%M/M/n/k server queue 
F5 = @(s,a,b,n,k) ...
    (b./(s+b))...
    .* (...
    sum(((a/b).^[0:n-1])./factorial([0:n-1])) ...
    + sum(((a/b).^[n:k-1]).*((n*b./(s+n*b)).^([1:k-n]))./(factorial(n).*n.^([0:k-n-1])))...
    )...
    / sum(((a/b).^[0:k-1])./(factorial(min(n,[0:k-1])).*(n.^max(0, [0:k-1]-n))));
pdf5 = @(x,a,b,n,k) matlabFunction(ilaplace(F5(s,a,b,n,k),s,x));
cdf5 = @(x,a,b,n,k) matlabFunction(ilaplace(F5(s,a,b,n,k)./s,s,x));

%Erlang queue
F6 = @(s,a,n) (F4(s,0,a)).^n;

F_RCH = F1; %RACH
F_NPU = F2; %User to Backhaul (Uplink network processor queue at the Base station)
F_NPD = F2; %Backhaul to User (Downlink network processor queue at the Base station)
F_BHU = F4; %Backhaul router (to the servers)
F_BHD = F4; %Backhaul router (from the servers)
F_TrU = F4; %Core router (to the Cloud)
F_TrD = F4; %Core router (from the Cloud)
F_MEC = F5; %MEC server
F_CLD = F5; %Cloud server

%=============================
%=============================
%=============================

lambda = lambda_in;

%iterative search of the Nash equilibrium point
 
alpha_bottom = 0;
alpha_top = 1;
alpha_delta = alpha_step; 
num_iter=0;
while  alpha_top - alpha_bottom > alpha_delta
    num_iter=num_iter+1;
    fprintf('Iteration #%d --- Checking interval [%g ,%g]\n', num_iter, alpha_bottom, alpha_top)
    alpha_mid = (alpha_bottom+alpha_top)/2;
    num_tests = 3*test_points; 
    for ac = 1 : num_tests
        alpha_ref = alpha_mid + (2 * ac - num_tests - 1) / (num_tests-1) *(alpha_top-alpha_bottom)/2;
        alpha_others = alpha_mid;
        alpha_eq = (alpha_ref + (num_users-1) * alpha_others) / num_users;
        policy_randalph_once_released()
        pfail_gand_mec(ac) = (1 - cdf_mec_at_timeout)*(1-ploss(lambda*alpha_eq/mu(8),nmec,kmec)) + ploss(lambda*alpha_eq/mu(8),nmec,kmec);
        pfail_gand_cloud(ac) = (1 - cdf_cloud_at_timeout)*(1-ploss(lambda*(1-alpha_eq)/mu(9),ncloud,kcloud)) + ploss(lambda*(1-alpha_eq)/mu(9),ncloud,kcloud);
        pfail_gand_ego(ac) = alpha_ref * pfail_gand_mec(ac) + (1-alpha_ref) * pfail_gand_cloud(ac);
    end
    pfail_test_dw = mean(pfail_gand_ego(1 : num_tests/3));
    pfail_test_sm = mean(pfail_gand_ego(num_tests/3+1 : 2*num_tests/3));
    pfail_test_up = mean(pfail_gand_ego(2*num_tests/3+1 : num_tests));
    fprintf('\n = = = = = = >> same: %.18f\t up: %.18f\t down: %.18f\n\n',pfail_test_sm,pfail_test_up,pfail_test_dw)
    if pfail_test_up < pfail_test_dw && pfail_test_up < pfail_test_sm
        alpha_bottom = alpha_mid;
    elseif pfail_test_dw < pfail_test_up && pfail_test_dw < pfail_test_sm
        alpha_top = alpha_mid;
    else
        alpha_bottom = alpha_mid - (alpha_top-alpha_bottom)/4;
        alpha_top = alpha_mid + (alpha_top-alpha_bottom)/4;
    end
end
alpha_nash = (alpha_top + alpha_bottom) / 2;
pfail_nash = pfail_gand_ego((num_tests+1)/2);

alpha_ref = alpha_nash;
alpha_others = alpha_nash;
alpha_eq = alpha_nash;
policy_randalph_once_released()
pfail_rand_nash = (1 - cdf_latency_at_timeout)*(1-loss(lambda_in,alpha_eq)) + loss(lambda_in,alpha_eq);
pfail_rand_mec_nash = (1 - cdf_mec_at_timeout)*(1-ploss(lambda*alpha_eq/mu(8),nmec,kmec)) + ploss(lambda*alpha_eq/mu(8),nmec,kmec);
pfail_rand_cloud_nash = (1 - cdf_cloud_at_timeout)*(1-ploss(lambda*(1-alpha_eq)/mu(9),ncloud,kcloud)) + ploss(lambda*(1-alpha_eq)/mu(9),ncloud,kcloud);


%now search for the best configuration (alpha is the same for all)

for a_l = 1:length(alpha_in)
    alpha_ref = alpha_in(a_l);
    alpha_others = alpha_ref;
    alpha_eq = (alpha_ref  + (num_users-1) * alpha_others) / num_users;
    policy_randalph_once_released()
    pfail_rand(a_l) = (1 - cdf_latency_at_timeout) * (1-loss(lambda_in,alpha_eq)) + loss(lambda_in,alpha_eq);
end

[pf_min pf_pos] = min(pfail_rand);

fprintf('Pfail (Coordinated) : %.6g at alpha = %g\n', pf_min, alpha_in(pf_pos))
fprintf('Pfail (Nash)        : %.6g at alpha = %g\n', pfail_rand_nash, alpha_nash)

perf = [pf_min, alpha_in(pf_pos), pfail_rand_nash, alpha_nash];


