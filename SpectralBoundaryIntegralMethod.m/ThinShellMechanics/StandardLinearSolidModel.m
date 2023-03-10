function [nbrev, viscousStress_prev, epsilbrev_prev] = ...
                    StandardLinearSolidModel(nbrev, viscousStress_prev, ...
                                             epsilbrev, epsilbrev_prev, ...
                                             mu_Mem, Tau, DT)
%% Membrane viscoelasticity using standard linear solid model
A = (4*Tau - DT)/(4*Tau + DT);
B = 4*mu_Mem/(4*Tau + DT);
viscousStress(:,:,1) = A*viscousStress_prev(:,:,1) + ...
                       B*(epsilbrev(:,:,1) - epsilbrev_prev(:,:,1)); % viscousStress_thet_thet
viscousStress(:,:,2) = A*viscousStress_prev(:,:,2) + ...
                       B*(epsilbrev(:,:,3) - epsilbrev_prev(:,:,3)); % viscousStress_phi_thet
viscousStress(:,:,3) = A*viscousStress_prev(:,:,3) + ...
                       B*(epsilbrev(:,:,2) - epsilbrev_prev(:,:,2)); % viscousStress_thet_phi
viscousStress(:,:,4) = A*viscousStress_prev(:,:,4) + ...
                       B*(epsilbrev(:,:,4) - epsilbrev_prev(:,:,4)); % viscousStress_phi_phi
viscousStress_prev = viscousStress;
epsilbrev_prev = epsilbrev;
nbrev = nbrev + viscousStress;