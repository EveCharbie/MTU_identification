function [g0,g1] = muscletendonequation(a,fiberLength,tendonLengthening,muscleparameter,umtLength)
optimalFiberLength = muscleparameter(:,1) ;
phi0 = muscleparameter(:,2) ;
maximalIsometricForce = muscleparameter(:,3) ;
tendonSlackLength = muscleparameter(:,4) ;
%% Muscle Tendon Architecture Equations
%%%%%%%%%%%%%
pennationAngle =  asin(optimalFiberLength .*  sin(phi0) ./ fiberLength) ; % get pennation angle
muscleLength = fiberLength .* cos(pennationAngle);

tendonLength = tendonSlackLength + tendonLengthening; %umtLength  - muscleLength; %TODO: safety if  tendonLength < tendonSlackLenght

normalizedTendonLength = tendonLength ./ tendonSlackLength;
normalizedFiberLength = fiberLength ./ optimalFiberLength;


%% Tendon Force Equation
%%%%%%%%%%%%%
% Tendon force-length (S1)
kT = 35; c1 = 0.200; c2 = 0.995; c3 = 0.250; % tendon parameters

normalizedTendonForcePart1 = normalizedTendonLength * 0 ;
normalizedTendonForcePart2 = c1 .* exp(kT .* (normalizedTendonLength - c2)) - c3;
normalizedTendonForce = normalizedTendonForcePart2 ;
% normalizedTendonForce = normalizedTendonForcePart1 + if_else(normalizedTendonLength < 1, 0, normalizedTendonForcePart2); % if normalized length under 0 the force = 0

% normalizedTendonForce = c1 .* exp(kT .* (normalizedTendonLength - c2)) - c3;
tendonForce= normalizedTendonForce .* maximalIsometricForce ;

%% Muscle Force Equations
%%%%%%%%%%%%%
% Active muscle force-length (S2)
b11 = 0.814483478343008 ; b21 = 1.055033428970575 ; b31 = 0.162384573599574 ; b41 = 0.063303448465465 ; % first Gaussian coefficents
b12 = 0.433004984392647 ; b22 = 0.716775413397760; b32 = -0.029947116970696 ; b42 = 0.200356847296188 ; % second Gaussian coefficents
b13 = 0.100 ; b23 = 1.000 ; b33 = 0.5 * sqrt(0.5) ; b43 = 0.000 ; % third Gaussian coefficents

num3 = normalizedFiberLength - b23;
den3 = b33 + b43 * normalizedFiberLength;
FMtilde3 = b13 * exp(-0.5 * num3.^2 ./ den3.^2);
num1 = normalizedFiberLength - b21;
den1 = b31 + b41 * normalizedFiberLength;
FMtilde1 = b11 * exp(-0.5 * num1.^2 ./ den1.^2);
num2 = normalizedFiberLength - b22;
den2 = b32 + b42 * normalizedFiberLength;
FMtilde2 = b12 * exp(-0.5 * num2.^2 ./ den2.^2);

normalizedMuscleActiveForceLength = FMtilde1 + FMtilde2 + FMtilde3;

MuscleActiveForceLength = a' .* normalizedMuscleActiveForceLength .* maximalIsometricForce ;

% Passive muscle force-length (S3)

kpe = 4.0 ; e0 = 0.6 ;
normalizedMusclePassiveForcePart1 = normalizedFiberLength * 0 ;
normalizedMusclePassiveForcePart2 = (exp(((kpe .* (normalizedFiberLength - 1))./e0)) - 1)./ (exp(kpe) - 1) ;
normalizedMusclePassiveForce = normalizedMusclePassiveForcePart1 + if_else(normalizedFiberLength < 1, 0, normalizedMusclePassiveForcePart2); % if normalized length under 0 the force = 0

musclePassiveForce = normalizedMusclePassiveForce .* maximalIsometricForce ;

% Muscle force-velocity (S4)
% d1 -0.318
% d2 -8.149
% d3 -0.374
% d4 0.886
normalizedMuscleForceVelocity = 1 ; % vitesse = 0

%% Forces function
normalizedMuscleForce = a' .* normalizedMuscleActiveForceLength .* normalizedMuscleForceVelocity + normalizedMusclePassiveForce ;
muscleForce = normalizedMuscleForce .* maximalIsometricForce ;


%% constraint function 
g0 = tendonForce - (cos(pennationAngle) .* muscleForce);
g1 = umtLength - (muscleLength + tendonLength);
end
