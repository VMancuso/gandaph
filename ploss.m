function y = ploss(rho,n,k)
if rho == 0
    y = 0;
else
    r = [0:k];
    den = sum(rho.^r./( factorial(min(n, r)) .* n .^ max(0, r-n)) );
    y = rho.^k ./ (factorial(n).*n^(k-n).*den);
    %(rho.^r./(factorial(min(n, r)).*n.^max(0, r-n)))/den
end
end


