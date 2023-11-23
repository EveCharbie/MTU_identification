%% Musculo skeletical model (Opensim equivalent)P2
[ret, name] = system('hostname');
if strcmp(name(1:9), '942-27984')
    addpath('C:\Users\Stage\Desktop\Neuromusculoskeletal Modeling\Casadi')
elseif strcmp(name(1:27), 'MacBook-Air-de-mickaelbegon')
    addpath('/Users/mickaelbegon/Downloads/casadi-3.6.3-osx64-matlab2018b/')
end
import casadi.*
                        %% Model definition 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "a" is neuromuscular activation (between 0 and 1)
% "q" is the spatial skeleton configuration (x, y, z, theta_hip, theta_knee, theta_ankle)
% "knownparameter" is subject segment length and muscle insertion position
% "musculoskeletal_states" = [q, known_parameter]
% "neuromusculoskeletal_state" = [a, q, known_parameter]
% "muscleTendonParameters" : musculoskeletal parameters that are generally assumed not to change: maximal isometric muscle force (Fom), optimal fiber length (ℓom), tendon slack length (ℓst), and pennation angle at optimal fiber length (φo).
% "rootedvariables" = (tendonLengthening, fiberLength)
% "all_state" = [neuromusculoskeletal_state, rootedvariables)


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

q = vertcat(x,y,z,theta_hip, theta_knee, theta_ankle ); % skeleton configuration

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

R_0_leg = R_0_thigh * R_thigh_leg; 
R_0_talus = R_0_leg * R_leg_talus; 
R_0_calc = R_0_talus * R_talus_calc; 




    %% Joint center
    %%%%%%%%%%%%%

% HJC = R_0_thigh(1:3, 4); % Hip
% KJC = R_0_leg(1:3, 4); % Knee
% AJC = R_0_talus(1:3, 4); % Knee
% TJC = Rototranslate(R_0_calc, clac) ; % toe
% CALC = R_0_calc(1:3, 4);  %calcaneum

% Muscle insertions in local frame
Local_Insertion_tibialis_anterior = SX.sym('Local_Insertion_tibialis_anterior',3);
Local_Origin_tibialis_anterior = SX.sym('Local_Origin_tibialis_anterior',3);
Local_ViaPoint_tibialis_anterior = SX.sym('Local_ViaPoint_tibialis_anterior',3);

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
    Local_Insertion_gastrocnemius ); %TODO: add viapoint

known_parameters = vertcat(segment_geometry, muscle_insertion, Local_ViaPoint_tibialis_anterior);

    %% Muscle origins and insertion in R0
    %%%%%%%%%%%%%

% Origin_tibialis_anterior = Rototranslate(R_0_leg, Local_Origin_tibialis_anterior);  
% Insertion_tibialis_anterior = Rototranslate(R_0_calc, Local_Insertion_tibialis_anterior);
% ViaPoint_tibialis_anterior = Rototranslate(R_0_leg, Local_ViaPoint_tibialis_anterior);
% 
% Origin_soleus = Rototranslate(R_0_leg,  Local_Origin_soleus) ;
% Insertion_soleus = Rototranslate(R_0_calc, Local_Insertion_soleus);
% 
% Origin_gastrocnemius = Rototranslate(R_0_thigh, Local_Origin_gastrocnemius) ;
% Insertion_gastrocnemius = Rototranslate(R_0_calc, Local_Insertion_gastrocnemius) ;
% 
% Origin = horzcat(Origin_tibialis_anterior, Origin_soleus, Origin_gastrocnemius);
% Insertion = horzcat(Insertion_tibialis_anterior, Insertion_soleus, Insertion_gastrocnemius);
% ViaPoint = horzcat(ViaPoint_tibialis_anterior) ; 
% Markers = horzcat(HJC, KJC, AJC, TJC,CALC); 

    %% Muscle origins and insertion in thigh
    %%%%%%%%%%%%%

    % rototranslation matrix 
% R_thigh_leg = Rototranslation_Rz(thigh(1),thigh(2),thigh(3), -theta_knee);
% R_leg_talus = Rototranslation_Rz(leg(1),leg(2),leg(3), -theta_ankle );
% R_talus_calc = Rototranslation_Rz(talus(1),talus(2),talus(3), 0 );
R_thigh_talus = R_thigh_leg * R_leg_talus; 
R_thigh_calc = R_thigh_talus * R_talus_calc; 

    % joint center 
HJC = [0; 0; 0]; % Hip
KJC = R_thigh_leg(1:3, 4); % Knee
AJC = R_thigh_talus(1:3, 4); % Ankle
TJC = Rototranslate(R_thigh_calc, clac) ; % toe
CALC = R_thigh_calc(1:3, 4);  %calcaneum

    % insertion origine muscle 
Origin_tibialis_anterior = Rototranslate(R_thigh_leg, Local_Origin_tibialis_anterior);  
Insertion_tibialis_anterior = Rototranslate(R_thigh_calc, Local_Insertion_tibialis_anterior);
ViaPoint_tibialis_anterior = Rototranslate(R_thigh_leg, Local_ViaPoint_tibialis_anterior);

Origin_soleus = Rototranslate(R_thigh_leg,  Local_Origin_soleus) ;
Insertion_soleus = Rototranslate(R_thigh_calc, Local_Insertion_soleus);

Origin_gastrocnemius = Local_Origin_gastrocnemius;
Insertion_gastrocnemius = Rototranslate(R_thigh_calc, Local_Insertion_gastrocnemius) ;

Origin = horzcat(Origin_tibialis_anterior, Origin_soleus, Origin_gastrocnemius);
Insertion = horzcat(Insertion_tibialis_anterior, Insertion_soleus, Insertion_gastrocnemius);
ViaPoint = horzcat(ViaPoint_tibialis_anterior) ; 
Markers = horzcat(HJC, KJC, AJC, TJC,CALC); 




%% 3. Functions about model geometry/kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
musculoskeletal_states = vertcat(q, known_parameters) ; % all parameters of musculoskeletal model 

ForwardKinematics = Function('ForwardKinematics', ...
    {musculoskeletal_states}, {Origin, Insertion, ViaPoint, Markers}, ...
    {'musculoskeletal_states'}, {'Origin', 'Insertion', 'ViaPoint', 'Markers'}) ;

% umtLength = sqrt(sum((Insertion - Origin).^2))';
% umtLength = vertcat(norm(Insertion_tibialis_anterior - ViaPoint_tibialis_anterior, 1) + ...
%                     norm(ViaPoint_tibialis_anterior - Origin_tibialis_anterior, 1), .... % Tibialis
%                     norm(Insertion_soleus - Origin_soleus, 1),... % Soleus
%                     norm(Insertion_gastrocnemius - Origin_gastrocnemius, 1)) ; % Gast 

umtLength = vertcat( sqrt(sum((Insertion_tibialis_anterior - ViaPoint_tibialis_anterior).^2)) + ...
                     sqrt(sum((ViaPoint_tibialis_anterior - Origin_tibialis_anterior).^2)) , ... % Tibialis
                     sqrt(sum((Insertion_soleus - Origin_soleus).^2)) , ... % Soleus
                     sqrt(sum((Insertion_gastrocnemius - Origin_gastrocnemius).^2)) ) ; % Gast

momentArm = jacobian(umtLength, q);

getUMTLength = Function('UMTLength', ...
    {musculoskeletal_states}, {umtLength}, ...
    {'musculoskeletal_states'}, {'umt_length'});
getMomentArm = Function('MomentArm', ...
    {musculoskeletal_states}, {momentArm}, ...
    {'musculoskeletal_states'}, {'moment_arm'});

tibiaLength = vertcat(sqrt(sum((KJC - AJC).^2))) ; % Gast 
getTIBIALength = Function('getTIBIALength', ...
    {musculoskeletal_states}, {tibiaLength}, ...
    {'musculoskeletal_states'}, {'umt_length'});
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

rootedvariables = vertcat(fiberLength, tendonLengthening);


    %% Muscle Tendon Architecture Equations
    %%%%%%%%%%%%%
% muscle tendon architectural parameter
% tendon length --> input in function of tendon force 
% fiber length --> input in funcction of muscle force 
% muscle length --> input in muscle tendon equilibrium => length tnendon =
% length muscle = length fiber * cos(pennation angle) 

tendonLength = tendonSlackLength + tendonLengthening; %umtLength  - muscleLength; %TODO: safety if  tendonLength < tendonSlackLenght
muscleLength = umtLength - tendonLength; %fiberLength .* cos(pennationAngle); 

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

MuscleActiveForceLength = a .* normalizedMuscleActiveForceLength .* maximalIsometricForce ; 

% Passive muscle force-length (S3)

kpe = 4.0 ; e0 = 0.6 ;
% normalizedMusclePassiveForcePart1 = normalizedFiberLength * 0 ;
normalizedMusclePassiveForcePart1 =  0 ;
normalizedMusclePassiveForcePart2 = (exp(((kpe .* (normalizedFiberLength - 1))./e0)) - 1)./ (exp(kpe) - 1) ;
% normalizedMusclePassiveForce = normalizedMusclePassiveForcePart1 + if_else(normalizedFiberLength < 1, 0, normalizedMusclePassiveForcePart2); % if normalized length under 0 the force = 0 
normalizedMusclePassiveForce = if_else(normalizedFiberLength < 1, ...
    normalizedMusclePassiveForcePart1, ...
    normalizedMusclePassiveForcePart2); % if normalized length under 0 the force = 0 
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
neuromusculoskeletal_state = vertcat(a, q, known_parameters) ; 

all_states = vertcat(neuromusculoskeletal_state,rootedvariables);

getTendonForce = Function('getTendonForce', ...
    {all_states, muscleTendonParameters}, {tendonForce}, ...
    {'all_states','muscle_tendon_parameters'}, {'tendon_force'}) ;

getMuscleForce = Function('getMuscleForce', ...
    {all_states, muscleTendonParameters}, {muscleForce}, ...
    {'all_states','muscle_tendon_parameters'}, {'muscle_force'}) ;

% getPennationAngle = Function('getPennationAngle', ...
%     {all_states, muscleTendonParameters}, {pennationAngle}, ...
%     {'all_states','muscle_tendon_parameters'}, {'pennation_angle'});

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


%% Computing Joint Moments and Angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Moment articulaire Mj(θ, t) = ∑ i=1 m (ri(θ) ⋅ Fimt(θ, t))
jointMoment = momentArm' * tendonForce ;
jointMoment = jointMoment(end,end) ;
getJointMoment = Function('getJointMoment', ...
    {all_states, muscleTendonParameters}, {jointMoment}, ...
    {'all_states', 'muscle_tendon_parameters'}, {'jointMoment'});

%% Muscle-tendon equilibrium 
p = vertcat(neuromusculoskeletal_state , muscleTendonParameters) ;
FT = SX.sym('Tendon_force', nMuscles);
FM = SX.sym('Muscle_force', nMuscles);
pennationAngle = SX.sym('Pennation_angle', nMuscles);

% constraint 
g3 = FT - tendonForce ; 
g4 = FM - muscleForce ; 
g5 = umtLength - (cos(pennationAngle) .* fiberLength + tendonLength) ; 
g6 = (optimalFiberLength .* sin(phi0)) - (fiberLength .* sin(pennationAngle)) ; 
g7 = FM .* cos(pennationAngle) - FT ; 

unknown  = vertcat(FT, FM, tendonLengthening, fiberLength, pennationAngle);
equilibriumError = Function('equilibriumError', {unknown, p}, {vertcat(g3, g4, g5, g6, g7)},{'x', 'p'}, {'residuals'}); 

opts_kinsol = struct("constraints", ones(15, 1), 'print_level',1); % 1 means >= 0 
opts_newton = struct("constraints", ones(15, 1), 'print_iteration',1); % 1 means >= 0 

equilibrateMuscleTendonN = rootfinder('equilibrateMuscleTendonN','newton',equilibriumError, opts_newton) ;
equilibrateMuscleTendon = rootfinder('equilibrateMuscleTendon','kinsol',equilibriumError, opts_kinsol) ;

%% Muscle-tendon equilibrium for all muscle 
% input 
LUMT = SX.sym('UMT_length', nMuscles);
SX.sym('Muscle_Activation', nMuscles);

% unknown 
FT = SX.sym('Tendon_force', nMuscles);
FM = SX.sym('Muscle_force', nMuscles);
pennationAngle = SX.sym('Pennation_angle', nMuscles); 
fiberLength ; 
tendonLength ; 

% constraint 
g3 = FT - tendonForce ; 
g4 = FM - muscleForce ; 
g5 = LUMT - (cos(pennationAngle) .* fiberLength + tendonLength) ; 
g6 = (optimalFiberLength .* sin(phi0)) - (fiberLength .* sin(pennationAngle)) ; 
g7 = FM .* cos(pennationAngle) - FT ; 

unknown  = vertcat(FT, FM, tendonLengthening, fiberLength, pennationAngle) ; 
known = vertcat(a, LUMT, muscleTendonParameters) ; 
equilibriumError1 = Function('equilibriumError1', {unknown, known}, {vertcat(g3, g4, g5, g6, g7)},{'x', 'p'}, {'residuals'}); 
equilibrateMuscleTendon1 = rootfinder('equilibrateMuscleTendon1','kinsol',equilibriumError1, opts_kinsol) ;

%% Muscle-tendon equilibrium for a single muscle 
unknown  = vertcat(FT(1), FM(1), tendonLengthening(1), fiberLength(1), pennationAngle(1)) ; 
known = vertcat(a(1), LUMT(1), optimalFiberLength(1), phi0(1),maximalIsometricForce(1),tendonSlackLength(1)) ; 
equilibriumErrorSingleMuscle = Function('equilibriumErrorSingleMuscle', {unknown, known}, {vertcat(g3(1), g4(1), g5(1), g6(1), g7(1))},{'x', 'p'}, {'residuals'}); 
opts_kinsol = struct("constraints", ones(5, 1), 'print_level',1); % 1 means >= 0 
equilibrateMuscleTendonSingleMuscle = rootfinder('equilibrateMuscleTendonSingleMuscle','kinsol',equilibriumErrorSingleMuscle, opts_kinsol) ;


%% representation function 
representationMusclePassiveForce = Function('representationMusclePassiveForce', ...
    {fiberLength(1), optimalFiberLength(1),maximalIsometricForce(1)}, {musclePassiveForce(1)}, ...
    {'fiberLength', 'optimalFiberLength','maximalIsometricForce'}, {'MusclePassiveForce'}) ;

representationMuscleActiveForceLength = Function('representationMuscleActiveForceLength', ...
    {a(1),fiberLength(1), optimalFiberLength(1),maximalIsometricForce(1)}, {MuscleActiveForceLength(1)}, ...
    {'a','fiberLength', 'optimalFiberLength','maximalIsometricForce'}, {'MuscleActiveForceLength'}) ;

representationTendonForce = Function('representationTendonForce', ...
    {tendonSlackLength(1), tendonLengthening(1),maximalIsometricForce(1)}, {tendonForce(1)}, ...
    {'tendonSlackLength', 'tendonLengthening','maximalIsometricForce'}, {'tendonForce'}) ;


%% struct that stock casadi function 
casadiFun = struct( ...
    'equilibrateMuscleTendon', equilibrateMuscleTendon,...
    'equilibrateMuscleTendon1', equilibrateMuscleTendon1,...
    'equilibrateMuscleTendonSingleMuscle', equilibrateMuscleTendonSingleMuscle,...
    'ForwardKinematics',ForwardKinematics,...
    'getUMTLength', getUMTLength,...
    'getMomentArm', getMomentArm,...
    'getJointMoment',getJointMoment,...
    'getTendonForce', getTendonForce, ...
    'getMuscleForce', getMuscleForce,...
    'getMusclePassiveForce',getMusclePassiveForce,...
    'getMuscleActiveForce', getMuscleActiveForce,...
    'normalizeTendonForce', normalizeTendonForce,...
    'normalizeTendonLength',normalizeTendonLength,...
    'normalizeFiberLength', normalizeFiberLength, ...
    'representationMusclePassiveForce', representationMusclePassiveForce, ...
    'representationMuscleActiveForceLength', representationMuscleActiveForceLength, ...
    'representationTendonForce', representationTendonForce);
 
