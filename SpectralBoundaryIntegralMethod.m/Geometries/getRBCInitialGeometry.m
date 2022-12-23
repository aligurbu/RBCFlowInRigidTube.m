function Xi = getRBCInitialGeometry(thet, phi, InitXi, InitOrient)
%%
% Undeformed shape of RBC taken from 
% Pozrikidis, C. (2005). "Axisymmetric motion of a file of red blood cells 
% through capillaries." Physics of Fluids, 17(3) [Equation (31)].
% (InitXi,InitOrient) give the initial position and orientation of the RBC
%%
nlat = size(thet,1);
nlon = size(phi,2);

alpha = 1.3858189;

Xi = zeros(nlat, nlon, 3);
Xi(:,:,1) = alpha*(sin(thet).*cos(phi));
Xi(:,:,2) = alpha*(sin(thet).*sin(phi));
Xi(:,:,3) = (0.5*alpha)*...
            (0.207+2.003*sin(thet).^2-1.123*sin(thet).^4).*cos(thet).* ...
                                                              ones(1,nlon);

if nargin > 2
    for m = 1:nlat
        for n = 1:nlon
            x = reshape(Xi(m,n,:),3,1);
            x = InitXi + InitOrient*x;
            Xi(m,n,:) = x;
        end
    end
end
