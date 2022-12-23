function mbrev = bendingStress(EB, kappabrev)
%% Bending stresses
mbrev(:,:,1) = EB*kappabrev(:,:,1);  % mbrev_thet_thet
mbrev(:,:,2) = EB*kappabrev(:,:,3);  % mbrev_phi_thet
mbrev(:,:,3) = EB*kappabrev(:,:,2);  % mbrev_thet_phi
mbrev(:,:,4) = EB*kappabrev(:,:,4);  % mbrev_phi_phi