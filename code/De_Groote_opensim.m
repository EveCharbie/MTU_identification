%% Musculo skeletical model (Opensim equivalent)

try 
addpath('/Users/mickaelbegon/Downloads/casadi-3.6.3-osx64-matlab2018b/')
catch 

end
try
 addpath('C:\Users\Stage\Desktop\Neuromusculoskeletal Modeling\Casadi')
catch

end 

import casadi.*



                        %% Musculoskeletal Geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Topology 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Definition 
    %%%%%%%%%%%%%
% B : Bones
% J : Joins 
% df : degree of freedom
     
% Total df = 6 df
% --> J : Hip (4 df : Rotx= 0 ; Roty = 0 ; Rotz = theta_hip ; X = x ; Y = y; Z = z )
%   --> B : Thigh
%    --> J : knee (1 df : Rotx= 0 ; Roty = 0 ; Rotz = theta_knee ; X = 0 ; Y = 0 ; Z = 0 )
%       --> B : Leg 
%           --> J : Ankle (1 dF: Rotx= 0 ; Roty = 0 ; Rotz = theta_ankle ; X = 0 ; Y = 0 ; Z = 0 )
%                --> B : Talus 
%                   --> J : Subtalar (0 dF: Rotx= 0 ; Roty = 0 ; Rotz = 0 ; X = 0 ; Y = 0 ; Z = 0 )
%                        --> B : Calc 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nMuscles = 3 ; NameMuscles = {'TibialisAnterior','Soleus','Gastrocnemius'} ; 
nBones = 4 ; NameBones = {'Thigh','Leg','Talus','Calcaneus'} ; 
nJoints = 4 ; NameJoints = {'Hip','Knee','Ankle','Subtalar'} ; 

%% 2. Casadi equation of the musculoskeletal model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% States (q's) 
    %%%%%%%%%%%%%

x = SX.sym('x');
y = SX.sym('y');
z = SX.sym('z');
theta_hip = SX.sym('theta_hip');
theta_knee = SX.sym('theta_knee');
theta_ankle = SX.sym('theta_ankle');

q = vertcat(x,y,z,theta_hip, theta_knee, theta_ankle );

    %% Model parameters that are personalized by scaling
    %%%%%%%%%%%%%

thigh = SX.sym('thigh',3);
leg = SX.sym('leg',3);
talus = SX.sym('talus',3);
clac = SX.sym('foot',3);

segment_geometry = vertcat(thigh, leg, talus ,clac);

    %% Rototanslation matrix
    %%%%%%%%%%%%%

R_0_thigh = Rototranslation_Rz(x,y,z, theta_hip);
R_thigh_leg = Rototranslation_Rz(thigh(1),thigh(2),thigh(3), -theta_knee);
R_leg_talus = Rototranslation_Rz(leg(1),leg(2),leg(3), -theta_ankle );
R_talus_calc = Rototranslation_Rz(talus(1),talus(2),talus(3), 0 );

R_0_leg = R_0_thigh * R_thigh_leg ; 
R_0_talus = R_0_leg * R_leg_talus ; 
R_0_calc = R_0_talus * R_talus_calc ; 

    %% Joint center
    %%%%%%%%%%%%%

HJC = R_0_thigh(1:3, 4); % Hip
KJC = R_0_leg(1:3, 4); % Knee
AJC = R_0_talus(1:3, 4); % Knee
TJC = Rototranslate(R_0_calc, clac) ; % toe
CALC = R_0_calc(1:3, 4);  %calcaneum

% Muscle insertions in local frame
Local_Insertion_tibialis_anterior = SX.sym('Local_Insertion_tibialis_anterior',3);
Local_Origin_tibialis_anterior = SX.sym('Local_Origin_tibialis_anterior',3);
Local_ViaPoint_tibialis_anterior = SX.sym('Local_Insertion_tibialis_anterior',3);

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

known_parameters = vertcat(segment_geometry, muscle_insertion,Local_ViaPoint_tibialis_anterior);

    %% Muscle origins and insertion in R0
    %%%%%%%%%%%%%

Origin_tibialis_anterior = Rototranslate(R_0_leg, Local_Origin_tibialis_anterior);  
Insertion_tibialis_anterior = Rototranslate(R_0_calc, Local_Insertion_tibialis_anterior);
ViaPoint_tibialis_anterior = Rototranslate(R_0_leg, Local_ViaPoint_tibialis_anterior);

Origin_soleus = Rototranslate(R_0_leg,  Local_Origin_soleus) ;
Insertion_soleus = Rototranslate(R_0_calc, Local_Insertion_soleus);
Origin_gastrocnemius = Rototranslate(R_0_thigh, Local_Origin_gastrocnemius) ;
Insertion_gastrocnemius = Rototranslate(R_0_calc, Local_Insertion_gastrocnemius) ;

Origin = horzcat(Origin_tibialis_anterior, Origin_soleus, Origin_gastrocnemius);
Insertion = horzcat(Insertion_tibialis_anterior, Insertion_soleus, Insertion_gastrocnemius);
ViaPoint = horzcat(ViaPoint_tibialis_anterior) ; 
Markers = horzcat(HJC, KJC, AJC, TJC,CALC); 

%% 3. Functions about model geometry/kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
musculoskeletal_states = vertcat(q, known_parameters) ; % all parameters of musculoskeletal model 

ForwardKinematics = Function('ForwardKinematics', ...
    {q, known_parameters}, {Origin, Insertion,ViaPoint, Markers}, ...
    {'q', 'known_parameters'}, {'Origin', 'Insertion','ViaPoint', 'Markers'}) ;

% umtLength = sqrt(sum((Insertion - Origin).^2))';
umtLength = vertcat( sqrt(sum((Insertion(:,1) - ViaPoint(:,1)).^2)) + sqrt(sum((ViaPoint(:,1) - Origin(:,1)).^2)),... % Tibialis
    sqrt(sum((Insertion(:,2) - Origin(:,2)).^2)),... % Soleus
    sqrt(sum((Insertion(:,3) - Origin(:,3)).^2))) ; % Gast 

momentArm = jacobian(umtLength, q);

getUMTLength = Function('UMTLength', ...
    {q, known_parameters}, {umtLength}, ...
    {'q', 'known_parameters'}, {'umt_length'});
getMomentArm = Function('MomentArm', ...
    {q, known_parameters}, {momentArm}, ...
    {'q', 'known_parameters'}, {'moment_arm'});



                        %% Muscle-tendon equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Activation Dynamics (Not use yet) 
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

% taua = 0.015; %SX.sym('Activation_time_constant',nMuscles);
% taud = 0.060; %SX.sym('Desactivation_time_constant',nMuscles);
% b = 0.1; %SX.sym('Transition_smoothness', nMuscles);
% 
% fa = 0.5 .* tanh(b .* (e-a)) ;
% da_dt = ((1 ./ taua .* (0.5 + 1.5 .* a)) .* (fa + 0.5) + ...
%     ((0.5 + 1.5 .* a) ./ (taud)) .* (-fa + 0.5)) .* (e - a);
% ode = struct('x',a, 'u', e, 'ode', da_dt);
% tgrid = linspace(0, 1, 100);
% activation_dynamics = integrator('activation_dynamics', 'cvodes', ode, tgrid(1), tgrid(1:end));
% % example: 
% next_a = activation_dynamics('x0',[0;0;0],'u',[0.5; 0.7; 1]);

%% 1. Muscle Contraction Dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Groote F De, Kinney AL, Rao A V, Fregly BJ. Evaluation of Direct Collocation
% Optimal Control Problem Formulations for Solving the Muscle Redundancy Problem.
% Ann Biomed Eng. 2016; 44: 2922–2936. https://doi.org/10.1007/s10439-016-1591-9
% PMID: 27001399

    %% Muscle Tendon Parameters (Fom, ℓom, ℓst, φo)
    %%%%%%%%%%%%%

% muscle tendon properties to be identified
optimalFiberLength = SX.sym('Optimal_fiber_length', nMuscles);
phi0 = SX.sym('Pennation_angle_at_muscle_optimal_fiber_length', nMuscles);
maximalIsometricForce = SX.sym('Maximal_isometric_muscle_force', nMuscles);
tendonSlackLength = SX.sym('Tendon_slack_length', nMuscles);

muscleTendonParameters =  vertcat(optimalFiberLength, phi0, maximalIsometricForce, tendonSlackLength) ; % TODO: should include more parameters? [MB]


% muscle tendon states input of the muscle tendon equation 
fiberLength = SX.sym('Fiber_length', nMuscles);
tendonLengthening = SX.sym('Tendon_Lengthening', nMuscles);
a = SX.sym('Muscle_Activation', nMuscles);

muscle_tendon_states = vertcat(fiberLength, tendonLengthening, a);

    %% Muscle Tendon Architecture Equations
    %%%%%%%%%%%%%
% muscle tendon architectural parameter
% tendon length --> input in function of tendon force 
% fiber length --> input in funcction of muscle force 
% muscle length --> input in muscle tendon equilibrium => length tnendon =
% length muscle = length fiber * cos(pennation angle) 


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
normalizedTendonForce = normalizedTendonForcePart1 + if_else(normalizedTendonLength < 1, 0, normalizedTendonForcePart2); % if normalized length under 0 the force = 0 

% normalizedTendonForce = c1 .* exp(kT .* (normalizedTendonLength - c2)) - c3; 
tendonForce= normalizedTendonForce .* maximalIsometricForce ; 


    %% Muscle Force Equations
    %%%%%%%%%%%%%
% Active muscle force-length (S2)
% b11 = 0.815 ; b21 = 1.055 ; b31 = 0.162 ;  b41 = 0.063 ; % first Gaussian coefficents
% b12 = 0.433 ; b22 = 0.717 ; b32 = -0.030 ; b42 = 0.200 ; % second Gaussian coefficents
% b13 = 0.100 ; b23 = 1.000 ; b33 = 0.354 ;  b43 = 0.000 ; % third Gaussian coefficents
% 
% normalizedMuscleActiveForceLength = (b11 .* exp((-0.5.* (normalizedFiberLength - b21).^2)./ (b31 + b41 .* normalizedFiberLength))) + ...
%     (b12 .* exp((-0.5.* (normalizedFiberLength - b22).^2)./ (b32 + b42 .* normalizedFiberLength))) + ...
%     (b13 .* exp((-0.5.* (normalizedFiberLength - b23).^2)./ (b33 + b43 .* normalizedFiberLength))) ;

[normalizedMuscleActiveForceLength] = MuscleForceLengthCSI(normalizedFiberLength,nMuscles);

MuscleActiveForceLength = normalizedMuscleActiveForceLength .* maximalIsometricForce ; 


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
normalizedMuscleForce = a .* normalizedMuscleActiveForceLength .* normalizedMuscleForceVelocity + normalizedMusclePassiveForce ;   
muscleForce = normalizedMuscleForce .* maximalIsometricForce ; 

%% Create all muscle-tendon functions
all_states = vertcat(musculoskeletal_states, muscle_tendon_states);

getTendonForce = Function('getTendonForce', ...
    {all_states, muscleTendonParameters}, {tendonForce}, ...
    {'all_states','muscle_tendon_parameters'}, {'tendon_force'}) ;

getMuscleForce = Function('getMuscleForce', ...
    {all_states, muscleTendonParameters}, {muscleForce}, ...
    {'all_states','muscle_tendon_parameters'}, {'muscle_force'}) ;

getPennationAngle = Function('getPennationAngle', ...
    {all_states, muscleTendonParameters}, {pennationAngle}, ...
    {'all_states','muscle_tendon_parameters'}, {'pennation_angle'});

normalizeTendonForce = Function('normalizeTendonForce', ...
    {all_states, muscleTendonParameters}, {normalizedTendonForce}, ...
    {'all_states', 'muscleTendonParameters'}, {'NormalizedTendonForce'});

getMusclePassiveForce = Function('getMusclePassiveForce', ...
    {all_states, muscleTendonParameters}, {musclePassiveForce}, ...
    {'all_states','muscle_tendon_parameters'}, {'MusclePassiveForce'}) ;

getMuscleActiveForce = Function('getMuscleActiveForce', ...
    {all_states, muscleTendonParameters}, {MuscleActiveForceLength}, ...
    {'all_states','muscle_tendon_parameters'}, {'MuscleActiveForce'}) ;

normalizeTendonLength = Function('normalizeTendonLength', ...
    {all_states, muscleTendonParameters}, {normalizedTendonLength}, ...
    {'all_states','muscle_tendon_parameters'}, {'length_tendon'});

normalizeFiberLength = Function('normalizeFiberLength', ...
    {all_states, muscleTendonParameters}, {normalizedFiberLength}, ...
    {'all_states','muscle_tendon_parameters'}, {'normalizeFiberLength'});


% TODO: not sure these functions are relevant [MB] according to me, those
% are not [AM]
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
    {all_states, muscleTendonParameters}, {jointMoment}, ...
    {'all_states', 'muscle_tendon_parameters'}, {'jointMoment'});


%% Muscle-tendon equilibrium
% determine muscle length such that TendonForce - (cos(pennation_angle) .* MuscleForce) = 0

g0 = tendonForce - (cos(pennationAngle) .* muscleForce); 
g1 = umtLength - (muscleLength + tendonLength); 

% g = Function('g',[tendonLength, FiberLength],[g0, g1]) ; 
x = vertcat(tendonLengthening, fiberLength);
p = vertcat(a,musculoskeletal_states , muscleTendonParameters) ;

g = Function('g', {x, p}, {vertcat(g0, g1)},{'x', 'p'}, {'residuals'}); 

opts = struct("constraints", ones(6,1)); % 1 means >= 0 
equilibrateMuscleTendon = rootfinder('equilibrateMuscleTendon','newton',g, opts) ; 



%% function d'optimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

unknown_parameters = horzcat(optimalFiberLength,phi0,maximalIsometricForce,tendonSlackLength) ;


%% struct that stock casadi function 
casadiFun = struct('ForwardKinematics',ForwardKinematics,...
    'getUMTLength', getUMTLength,...
    'getMomentArm', getMomentArm,...
    'getTendonForce', getTendonForce, ...
    'getMuscleForce', getMuscleForce,...
    'getPennationAngle', getPennationAngle,...
    'getMusclePassiveForce',getMusclePassiveForce,...
    'getMuscleActiveForce', getMuscleActiveForce,...
    'normalizeTendonForce', normalizeTendonForce,...
    'normalizeTendonLength',normalizeTendonLength,...
    'normalizeFiberLength', normalizeFiberLength);
 
