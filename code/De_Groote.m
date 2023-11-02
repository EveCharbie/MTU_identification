%% Casadi
%addpath('/Users/mickaelbegon/Downloads/casadi-3.6.3-osx64-matlab2018b/')
addpath('C:\Users\Stage\Desktop\Neuromusculoskeletal Modeling\Casadi')

import casadi.*



%% Musculoskeletal Geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% definition du model
% 3 Segment (bones): foot, leg thigh
% 3 Muscles: Tibialis anterior, soleus, gastrocnemius

% Joint Center :
% TJC : toe joint center
% AJC : ankle joint center
% KJC : knee joint center
% HJC : hip joint center

%% States
x = SX.sym('x');
y = SX.sym('y');
z = SX.sym('z');
theta_foot = SX.sym('theta_foot');
theta_ankle = SX.sym('theta_ankle');
theta_knee = SX.sym('theta_knee');
q = vertcat(x,y,z,theta_foot, theta_ankle, theta_knee);

nMuscles = 3;
e = SX.sym('Neuromuscular_activation',nMuscles);
a = SX.sym('Muscle_activation',nMuscles);
muscleLength = SX.sym('muscle_length',nMuscles); %TODO: decide with length to use [MB] and adapt code accordingly

states = vertcat(x,y,z,theta_foot, theta_ankle, theta_knee, a, muscleLength);


%% Model parameters that are personalized by scaling
length_foot = SX.sym('length_foot');
length_leg = SX.sym('length_leg');
length_thigh = SX.sym('length_thigh');
segment_length = vertcat(length_foot, length_leg, length_thigh);

% Muscle insertions in local frame
Local_Insertion_tibialis_anterior = SX.sym('Local_Insertion_tibialis_anterior',3);
Local_Origin_tibialis_anterior = SX.sym('Local_Origin_tibialis_anterior',3);
Local_Insertion_soleus = SX.sym('Local_Insertion_soleus',3);
Local_Origin_soleus = SX.sym('Local_Origin_soleus',3);
Local_Insertion_gastrocnemius = SX.sym('Local_Insertion_gastrocnemius',3);
Local_Origin_gastrocnemius = SX.sym('Local_Origin_gastrocnemius',3);
muscle_insertion = vertcat( ...
    Local_Origin_tibialis_anterior, ...
    Local_Origin_soleus, ...
    Local_Origin_gastrocnemius, ...
    Local_Insertion_tibialis_anterior, ...
    Local_Insertion_soleus, ...
    Local_Insertion_gastrocnemius );

known_parameters = vertcat(segment_length, muscle_insertion);


%% model geometric
% Rototranslation
% longitudinal axis:  x

R_0_foot = Rototranslation_Rz(x,y,z, theta_foot);
R_foot_leg = Rototranslation_Rz(-length_foot,0,0, theta_ankle-pi/2);
R_leg_thigh = Rototranslation_Rz(-length_leg,0,0, theta_knee);

R_0_leg = R_0_foot * R_foot_leg; 
R_0_thigh = R_0_leg * R_leg_thigh;

%% Joint center
TJC = R_0_foot(1:3, 4);  % Toe
AJC = R_0_leg(1:3, 4);   % Ankle
KJC = R_0_thigh(1:3, 4); % Knee
HJC = Rototranslate(R_0_thigh,  [-length_thigh; 0; 0;]) ; % Hip

%% muscle origins and insertion in R0
Origin_tibialis_anterior = Rototranslate(R_0_leg, Local_Origin_tibialis_anterior);  
Insertion_tibialis_anterior = Rototranslate(R_0_foot, Local_Insertion_tibialis_anterior);
Origin_soleus = Rototranslate(R_0_leg,  Local_Origin_soleus) ;
Insertion_soleus = Rototranslate(R_0_foot, Local_Insertion_soleus);
Origin_gastrocnemius = Rototranslate(R_0_thigh, Local_Origin_gastrocnemius) ;
Insertion_gastrocnemius = Rototranslate(R_0_foot, Local_Insertion_gastrocnemius) ;

Origin = horzcat(Origin_tibialis_anterior, Origin_soleus, Origin_gastrocnemius);
Insertion = horzcat(Insertion_tibialis_anterior, Insertion_soleus, Insertion_gastrocnemius);
Markers = horzcat(TJC, AJC, KJC, HJC);


%% functions about model geometry/kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ForwardKinematics = Function('ForwardKinematics', ...
    {q, known_parameters}, {Origin, Insertion, Markers}, ...
    {'q', 'known_parameters'}, {'Origin', 'Insertion', 'Markers'}) ;

umtLength = sqrt(sum((Insertion - Origin).^2))';
momentArm = jacobian(umtLength, q);

getUMTLength = Function('UMTLength', ...
    {q, known_parameters}, {umtLength}, ...
    {'q', 'known_parameters'}, {'umt_length'});
getMomentArm = Function('MomentArm', ...
    {q, known_parameters}, {momentArm}, ...
    {'q', 'known_parameters'}, {'moment_arm'});



%% Activation Dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Groote F De, Kinney AL, Rao A V, Fregly BJ. Evaluation of Direct Collocation
% Optimal Control Problem Formulations for Solving the Muscle Redundancy Problem.
% Ann Biomed Eng. 2016; 44: 2922–2936. https://doi.org/10.1007/s10439-016-1591-9
% PMID: 27001399

% muscle activation is described by two nonlinear, first ODE

% parameter value
% taua = 0.015
% taud = 0.060
% b = 0.1

taua = SX.sym('Activation_time_constant',nMuscles);
taud = SX.sym('Desactivation_time_constant',nMuscles);
b = SX.sym('Transition_smoothness', nMuscles);

fa = 0.5 .* tanh(b .* (e-a)) ;
da_dt = ((1 ./ taua .* (0.5 + 1.5 .* a)) .* (fa + 0.5) + ...
    ((0.5 + 1.5 .* a) ./ (taud)) .* (-fa + 0.5)) .* (e - a);


%% Muscle Contraction Dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Groote F De, Kinney AL, Rao A V, Fregly BJ. Evaluation of Direct Collocation
% Optimal Control Problem Formulations for Solving the Muscle Redundancy Problem.
% Ann Biomed Eng. 2016; 44: 2922–2936. https://doi.org/10.1007/s10439-016-1591-9
% PMID: 27001399

% muscle tendon properties to be identified
optimalFiberLength = SX.sym('Optimal_fiber_length', nMuscles);
phi0 = SX.sym('Pennation_angle_at_muscle_optimal_fiber_length', nMuscles);
maximalIsometricForce = SX.sym('Maximal_isometric_muscle_force', nMuscles);
tendonSlackLength = SX.sym('Tendon_slack_length', nMuscles);

muscleTendonParameters =  vertcat(optimalFiberLength, phi0, maximalIsometricForce, tendonSlackLength) ; % TODO: should include more parameters? [MB]


%% architectural parameter [AM]
% muscle tendon architectural parameter
% tendon length --> input in function of tendon force 
% fiber length --> input in funcction of muscle force 
% muscle length --> input in muscle tendon equilibrium => length tnendon =
% length muscle = length fiber * cos(pennation angle) 


pennationAngle =  asin(optimalFiberLength .*  sin(phi0) ./ fiberLength) ; % get pennation angle
muscleLength = fiberLength * cos(pennationAngle) ; 
tendonLength = umtLength  - muscleLength; %TODO: safety if  tendonLength < tendonSlackLenght

normalizedTendonLength = tendonLength ./ tendonSlackLength; 
normalizedFiberLength = fiberLength ./ optimalFiberLength; 



%% Tendon
% Tendon force-length (S1)
kT = 35; c1 = 0.200; c2 = 0.995; c3 = 0.250; % tendon parameters
normalizedTendonForce = c1 .* exp(kT .* (normalizedTendonLength - c2)) - c3; 
tendonForce= normalizedTendonForce .* maximalIsometricForce ; 


%% Muscle
% Active muscle force-length (S2)
b11 = 0.815 ; b21 = 1.055 ; b31 = 0.162 ;  b41 = 0.063 ; % first Gaussian coefficents
b12 = 0.433 ; b22 = 0.717 ; b32 = -0.030 ; b42 = 0.200 ; % second Gaussian coefficents
b13 = 0.100 ; b23 = 1.000 ; b33 = 0.354 ;  b43 = 0.000 ; % third Gaussian coefficents

normalizedMuscleActiveForceLength = (b11 .* exp((-0.5.* (normalizedFiberLength - b21).^2)./ (b31 + b41 .* normalizedFiberLength))) + ...
    (b12 .* exp((-0.5.* (normalizedFiberLength - b22).^2)./ (b32 + b42 .* normalizedFiberLength))) + ...
    (b13 .* exp((-0.5.* (normalizedFiberLength - b23).^2)./ (b33 + b43 .* normalizedFiberLength))) ;

MuscleActiveForceLength = normalizedMuscleActiveForceLength .* maximalIsometricForce ; 

% Passive muscle force-length (S3)
kpe = 4.0 ; e0 = 0.6 ;
normalizedMusclePassiveForce = (exp(((kpe .* (normalizedFiberLength - 1))./e0)) - 1)./ (exp(kpe) - 1) ;


musclePassiveForce = normalizedMusclePassiveForce .* maximalIsometricForce ; 

% Muscle force-velocity (S4)
% d1 -0.318
% d2 -8.149
% d3 -0.374
% d4 0.886
NormalizedMuscleForceVelocity = 1 ; % vitesse = 0



%% Forces function
normalizedMuscleForce = a .* normalizedMuscleActiveForceLength .* NormalizedMuscleForceVelocity + normalizedMusclePassiveForce ;   
muscleForce = normalizedMuscleForce .* maximalIsometricForce ; 
1+1 

%% Create all muscle-tendon functions
getTendonForce = Function('getTendonForce', ...
    {states, known_parameters, muscleTendonParameters}, {tendonForce}, ...
    {'states','known_parameters','muscle_tendon_parameters'}, {'tendon_force'}) ;

getMuscleForce = Function('getMuscleForce', ...
    {states, known_parameters, muscleTendonParameters}, {muscleForce}, ...
    {'states','known_parameters','muscle_tendon_parameters'}, {'muscle_force'}) ;

getPennationAngle = Function('getPennationAngle', ...
    {states, muscleTendonParameters}, {pennationAngle}, ...
    {'states','muscle_tendon_parameters'}, {'pennation_angle'});

normalizeTendonForce = Function('normalizeTendonForce', ...
    {states, known_parameters, muscleTendonParameters}, {normalizedTendonForce}, ...
    {'states','known_parameters', 'muscleTendonParameters'}, {'NormalizedTendonForce'});

getMusclePassiveForce = Function('getMusclePassiveForce', ...
    {states, known_parameters, muscleTendonParameters}, {musclePassiveForce}, ...
    {'states','known_parameters','muscle_tendon_parameters'}, {'MusclePassiveForce'}) ;

getMuscleActiveForce = Function('getMuscleActiveForce', ...
    {states, known_parameters, muscleTendonParameters}, {MuscleActiveForceLength}, ...
    {'states','known_parameters','muscle_tendon_parameters'}, {'MuscleActiveForce'}) ;


% TODO: not sure these functions are relevant [MB]
% normalizeTendonLength = Function('normalizeTendonLength', ...
%     {tendonLength, muscleTendonParameters}, {normalizedTendonLength}, ...
%     {'length_tendon_normalized','muscle_tendon_parameters'}, {'length_tendon'});
% 
% normalizeFiberLength = Function('normalizeFiberLength', ...
%     {FiberLength, muscleTendonParameters}, {normalizedFiberLength}, ...
%     {'FiberLength','muscleTendonParameters'}, {'normalizeFiberLength'});
% 
% getNormalizeMusclePassiveForce = Function('Normalized_Muscle_Passive_Force_Length', ...
%     {states, known_parameters}, {normalizedMusclePassiveForce}, ...
%     {'states','known_parameters'}, {'NormalizedMusclePassiveForceLength'});
% 
% getNormalizedMuscleActiveForceLength = Function('Normalized_Muscle_Active_Force_Length', ...
%     {states,known_parameters}, {normalizedMuscleActiveForceLength}, ...
%     {'states','known_parameters'}, {'NormalizedMuscleActiveForceLength'}) ;
% 
% getNormalizedMuscleForce = Function('getNormalizedMuscleForce', ...
%     {states, known_parameters}, {normalizedMuscleForce}, ...
%     {'states', 'known_parameters'}, {'normalizedMuscleForce'});

%% Computing Joint Moments and Angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Moment articulaire Mj(θ, t) = ∑ i=1 m (ri(θ) ⋅ Fimt(θ, t))
jointMoment = momentArm' * tendonForce ;
getJointMoment = Function('getJointMoment', ...
    {states, known_parameters, muscleTendonParameters}, {jointMoment}, ...
    {'states', 'known_parameters', 'muscle_tendon_parameters'}, {'jointMoment'});


%% Muscle-tendon equilibrium
% determine muscle length such that TendonForce - (cos(pennation_angle) .* MuscleForce) = 0

g0 = tendonForce - (cos(pennationAngle) .* muscleForce); 
% g1 = umtLength' - (cos(pennationAngle) .* FiberLength + tendonLength) ; 
% g = Function('g',[tendonLength, FiberLength],[g0, g1]) ; 
g = Function('g', [muscleLength, vertcat(states, known_parameters, muscleTendonParameters)], [g0]); 

equilibrateMuscleTendon = rootfinder('equilibrateMuscleTendon','newton',g) ; 




%% function d'optimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

unknown_parameters = horzcat(optimalFiberLength,phi0,maximalIsometricForce,tendonSlackLength) ;
%%
% f minimisation quadratique
T = SX.sym('Mesured_torque',1) ;

f = (T -sum(jointMoment)).^2 ;
% g = Ft - Fm*cos(pennation_angle)

% %%
% nlp = struct('x',unknown_parameters, 'f',f, 'g',g);
% S = nlpsol('S', 'ipopt', nlp);
% sol = S('x0',x0, 'lbx',min_param,'ubx',max_param, 'lbg',0,'ubg',0); %il faut que lbg soit un vecteur de 0 de la taille de la contrainte
% Param = sol.x;
% disp(x_opt)