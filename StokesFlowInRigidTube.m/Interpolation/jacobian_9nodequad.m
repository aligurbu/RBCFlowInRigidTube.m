function [normal_vec, JJ] = jacobian_9nodequad(s1, s2, xi_e)
%%
% Compute surface jacobian = ||dxi/ds1 X dxi/ds2|| 
% for 9-node quadrilateral element 
% at an intrinsic point (s1, s2) \in [-1, 1] x [-1, 1]
% xi_e = 3x9 matrix of node coordinates of element
%%
DN = deriv_interpolate_9nodequad(s1, s2);

dxids = xi_e*DN;

% cross product - for efficiency
normal_vec = [dxids(2,1)*dxids(3,2) - dxids(3,1)*dxids(2,2); ...
              dxids(3,1)*dxids(1,2) - dxids(1,1)*dxids(3,2); ...
              dxids(1,1)*dxids(2,2) - dxids(2,1)*dxids(1,2)];

JJ = norm(normal_vec);      % Jacobian 
normal_vec = normal_vec/JJ; % Unit normal vector