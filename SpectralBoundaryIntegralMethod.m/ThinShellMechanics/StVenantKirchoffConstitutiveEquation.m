function nbrev = StVenantKirchoffConstitutiveEquation(ED, ES, epsilbrev)
%% St. Venant Kirchhoff model
nbrev(:,:,1) = (ED+ES)*epsilbrev(:,:,1) + ED*epsilbrev(:,:,4); % nbrev_thet_thet
nbrev(:,:,2) = ES*epsilbrev(:,:,3);                            % nbrev_phi_thet
nbrev(:,:,3) = ES*epsilbrev(:,:,2);                            % nbrev_thet_phi
nbrev(:,:,4) = (ED+ES)*epsilbrev(:,:,4) + ED*epsilbrev(:,:,1); % nbrev_phi_phi