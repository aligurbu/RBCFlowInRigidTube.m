%% Reference parameters 
CellVolume = 93.78*10^(-18); % m^3
RefLength = ((3*CellVolume)/(4*pi))^(1/3); % m
RefShearRate = 100; % 1/s
RefVelocity = RefLength*RefShearRate; % m/s
RefViscosity = 1.2*10^(-3); % Pa.s = N.s/m^2
RefPressure = RefViscosity*RefShearRate; % Pa = N/m^2
RefElasticModulus = RefViscosity*RefVelocity; %N/m
RefBendingModulus = RefLength^2*RefViscosity*RefVelocity; % N.m