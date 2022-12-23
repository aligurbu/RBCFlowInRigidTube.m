function N = interpolate_9nodequad(s1, s2)
%%
% Compute interpolation functions for 9-node quadrilateral element 
% at an intrinsic point (s1, s2) \in [-1, 1] x [-1, 1]
%%
N = [(s1*(s1 - 1)*s2*(s2 - 1))/4;
     (s1*(s1 + 1)*s2*(s2 - 1))/4;
     (s1*(s1 + 1)*s2*(s2 + 1))/4;
     (s1*(s1 - 1)*s2*(s2 + 1))/4;
     ((1 - s1^2)*s2*(s2 - 1))/2;
     (s1*(s1 + 1)*(1 - s2^2))/2;
     ((1 - s1^2)*s2*(s2 + 1))/2;
     (s1*(s1 - 1)*(1 - s2^2))/2;
     (1 - s1^2)*(1 - s2^2)];