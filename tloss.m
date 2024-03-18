function y = tloss(lambda,alp,mumec,nmec,kmec,mucloud,ncloud,kcloud)
y = alp*ploss(lambda*alp/mumec, nmec, kmec)+(1-alp)*ploss(lambda*(1-alp)/mucloud, ncloud, kcloud); 
end