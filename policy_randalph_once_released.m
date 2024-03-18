printing = 0; %do not print details

if alpha_ref == alpha_others
    fprintf('lambda = %g\talpha = %g\n',lambda_in,alpha_ref);
else
    fprintf('lambda = %g\talphas = %g %g %g\n',lambda_in,alpha_ref, alpha_others, alpha_eq);
end

%=======
%check the load of nodes with infinite waiting room
%
bf = zeros(1,9);
arr_rates = [lambda lambda lambda*(1-loss(lambda,alpha_eq))...
    lambda lambda*(1-loss(lambda,alpha_eq)) ...
    lambda*(1-alpha_eq) ...
    lambda*(1-alpha_eq)*(1-ploss(lambda*(1-alpha_eq)/mu(9),ncloud,kcloud))...
    lambda*alpha_eq ...
    lambda*(1-alpha_eq)
    ];
nodes = ["RACH   " "NetP-UL" "NetP-DL" "BaHa-UL" "BaHa-DL" "Core-UL" "Core-DL" "MEC    " "Cloud  "];

for j = 2:7
    if arr_rates(j)/mu(j) >= 1
        bf(j) = 1;
        if printing
            fprintf('%s: %9g / %9g <---\n', nodes(j), arr_rates(j), mu(j))
        end
    else
        if printing
            fprintf('%s: %9g / %9g\n', nodes(j), arr_rates(j), mu(j))
        end
    end
end
if printing
    fprintf('%s: %9g / %9g\n', nodes(8), arr_rates(8), mu(8)*nmec)
    fprintf('%s: %9g / %9g\n', nodes(9), arr_rates(9), mu(9)*ncloud)
    fprintf('loss: %g | mec: %g | cloud: %g\n',loss(lambda,alpha_eq),ploss(lambda*alpha_eq/mu(8),nmec,kmec),ploss(lambda*(1-alpha_eq)/mu(9),ncloud,kcloud))
    if bf(1) > 0
        fprintf('WARNING, THE RACH CANNOT BE LINEARIZED\n')
    end

    if sum(bf(2:end)) > 0
        fprintf('WARNING, UNSTABLE QUEUE(S)\n')
    end
end

%compute state probabilities of MEC and Cloud servers

pimec = ones(1,kmec+1);
rhomec = alpha_eq*lambda / mu(8);
for i = 2 : (nmec+1)
    pimec(i) = pimec(i-1) * rhomec / (i-1);
end
for i = (nmec + 2) : (kmec+1)
    pimec(i) = pimec(i-1) * rhomec / nmec;
end
pimec(1) = 1 / sum(pimec);
pimec(2:end) = pimec(2:end) * pimec(1);

picloud = ones(1,kcloud+1);
rhocloud = (1-alpha_eq)*lambda / mu(9);
for i = 2 : (ncloud+1)
    picloud(i) = picloud(i-1) * rhocloud / (i-1);
end
for i = (ncloud + 2) : (kcloud+1)
    picloud(i) = picloud(i-1) * rhocloud / ncloud;
end
picloud(1) = 1 / sum(picloud);
picloud(2:end) = picloud(2:end) * picloud(1);


rmec = lambda * alpha_eq / mu(8);
lmec = lambda * alpha_eq * (1-ploss(rmec,nmec,kmec));
rcloud = lambda * (1-alpha_eq) / mu(9);
lcloud = lambda * (1-alpha_eq) * (1-ploss(rcloud,ncloud,kcloud));

%compute the prob of timeout by inverting the LST of the latency 
cdf_mec_at_timeout = euler_inversion(matlabFunction(...
    F_RCH(s)...
    .*F_NPU(s,lambda,mu(2))...
    .*F_BHU(s,lambda,mu(4))...
    .*F_MEC(s,lambda*alpha_eq,mu(8),nmec,kmec)...
    .*F_BHD(s,lambda*(1-loss(lambda,alpha_eq)),mu(5))...
    .*F_NPD(s,lambda*(1-loss(lambda,alpha_eq)),mu(3))...
    ./s),timeout-del_mec);

cdf_cloud_at_timeout = euler_inversion(matlabFunction(...
    F_RCH(s)...
    .*F_NPU(s,lambda,mu(2))...
    .*F_BHU(s,lambda,mu(4))...
    .*F_TrU(s,lambda*(1-alpha_eq),mu(6))...
    .*F_CLD(s,lambda*(1-alpha_eq),mu(9),ncloud,kcloud)...
    .*F_TrD(s,lambda*(1-alpha_eq)...
    *(1-ploss(rcloud,ncloud,kcloud)),mu(7))...
    .*F_BHD(s,lambda*(1-loss(lambda,alpha_eq)),mu(5))...
    .*F_NPD(s,lambda*(1-loss(lambda,alpha_eq)),mu(3))...
    ./s),timeout-del_cloud);

cdf_latency_at_timeout = ...
    (lmec * cdf_mec_at_timeout + lcloud * cdf_cloud_at_timeout) / (lmec + lcloud);




