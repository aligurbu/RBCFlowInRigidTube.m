function [NGthet, NGphi, eta, wg] = setup_integration_grid(N, NGSphere)
%%
NGthet = NGSphere + 1;
NGphi = 2*NGSphere + 1;
[thetg,wg] = gaqdm(NGthet);

Psum = zeros(N+1,NGthet);
for n = 0:N
    P__ = legendre(n,cos(thetg));
    Psum(n+1,:) = P__(1,:);
end
Psum = cumsum(Psum,1);
Psum = Psum(end,:)';

eta = 2*sin(thetg/2).*Psum;