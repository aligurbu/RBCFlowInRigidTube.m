function zetamin = closest_zeta_9nodequad(chi, xi_e, dmin, dind)

r = chi - xi_e(:,dind);
[zeta_1, zeta_2] = getZeta(dind);
DN = deriv_interpolate_9nodequad(zeta_1, zeta_2);
dxidzeta = xi_e*DN;
zeta_n_1 = zeta_1 + (sum(r.*dxidzeta(:,1)))/norm(dxidzeta(:,1));
zeta_n_2 = zeta_2 + (sum(r.*dxidzeta(:,2)))/norm(dxidzeta(:,2));
if abs(zeta_n_1)>1
    zeta_n_1 = zeta_n_1/abs(zeta_n_1);
end
if abs(zeta_n_2)>1
    zeta_n_2 = zeta_n_2/abs(zeta_n_2);
end

N = interpolate_9nodequad(zeta_n_1, zeta_n_2);
r = chi - xi_e*N;
d_n = norm(r);

Error = abs(d_n - dmin)/abs(dmin);
while Error > 10^(-10)
    dmin = d_n;
    zeta_1 = zeta_n_1; zeta_2 = zeta_n_2; 
    DN = deriv_interpolate_9nodequad(zeta_1, zeta_2);
    dxidzeta = xi_e*DN;
    zeta_n_1 = zeta_1 + (sum(r.*dxidzeta(:,1)))/norm(dxidzeta(:,1));
    zeta_n_2 = zeta_2 + (sum(r.*dxidzeta(:,2)))/norm(dxidzeta(:,2));
    if abs(zeta_n_1)>1
        zeta_n_1 = zeta_n_1/abs(zeta_n_1);
    end
    if abs(zeta_n_2)>1
        zeta_n_2 = zeta_n_2/abs(zeta_n_2);
    end
    N = interpolate_9nodequad(zeta_n_1, zeta_n_2);
    r = chi - xi_e*N;
    d_n = norm(r);
    Error = abs(d_n - dmin)/abs(dmin);
end
zetamin = [zeta_1; zeta_2];
% dmin = d_n;
end
    