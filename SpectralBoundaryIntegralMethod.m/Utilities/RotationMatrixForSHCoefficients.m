function R_ = RotationMatrixForSHCoefficients(N, nlat, nlon, thet, phi)
%% Create rotation matrices
P = getSPcoeffTrans(N);
permut = [1:(N+1)*(N+2)/2 ...
         (N+1)^2+(1:(N+1)*(N+2)/2) ...
         2*(N+1)^2+(1:(N+1)*(N+2)/2) ...
         ((N+1)*(N+2)/2+1):(N+1)^2 ...
         (N+1)^2+(((N+1)*(N+2)/2+1):(N+1)^2) ...
         2*(N+1)^2+(((N+1)*(N+2)/2+1):(N+1)^2)];
R_ = cell(nlat,nlon);
indR = cumsum([1 (2*(0:N)+1).^2]);
nnzR = indR(end)-1; % number of non-zeros elements in R
IR = zeros(nnzR,1); % row indices of non-zero elements in R
JR = zeros(nnzR,1); % column indices of non-zero elements in R
for n = 0:N
    ind_ = repmat(((n^2+1):(n+1)^2)',1,2*n+1);
    indT_ = ind_';
    IR(indR(n+1):(indR(n+2)-1)) = ind_(:);
    JR(indR(n+1):(indR(n+2)-1)) = indT_(:);
end
for m = 1:nlat
    thetChi_ = thet(m);
    for n = 1:nlon
        phiChi_ = phi(n);
        r = expm( hat([-sin(phiChi_) cos(phiChi_) 0]*thetChi_) );
        SHR = SHRotate(r, N);
        SHRot = zeros(nnzR,1);
        for k = 0:N
            SHRot(indR(k+1):(indR(k+2)-1)) = SHR{k+1}(:);
        end
        Rsparse = sparse(IR, JR, SHRot);
        T = P\(Rsparse*P);
        Texpanded = blkdiag(T,T,T);
        R_{m,n} = Texpanded(permut,permut);
    end
end