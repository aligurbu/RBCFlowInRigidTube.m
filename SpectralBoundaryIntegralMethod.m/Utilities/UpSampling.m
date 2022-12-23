function cup = UpSampling(c,N,Nup)
%% Up-sample 

cup = zeros(Nup+1, Nup+1, size(c,3)); 
cup(1:N+1,1:N+1,:) = c;