%% Casadi
addpath('/Users/mickaelbegon/Downloads/casadi-3.6.3-osx64-matlab2018b/')
import casadi.*



%% Musculoskeletal Geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% definition du model
% Segment (bones): foot, leg thigh
% Muscles 

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

e = SX.sym('Neuromuscular_activation',3);
a = SX.sym('Muscle_activation',3);

states = vertcat(x,y,z,theta_foot, theta_ankle, theta_knee, a);


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
muscle_insertion = vertcat( Local_Origin_tibialis_anterior, ...
    Local_Origin_soleus, ...
    Local_Origin_gastrocnemius, ...
    Local_Insertion_tibialis_anterior, ...
    Local_Insertion_soleus, ...
    Local_Insertion_gastrocnemius );

known_parameters = vertcat(segment_length, muscle_insertion) ;


%% model geometric
% Rototranslation
% axe longitudinal x

R_0_foot = Rototranslation_Rz(x,y,z, theta_foot);
R_foot_leg = Rototranslation_Rz(-length_foot,0,0, theta_ankle-pi/2);
R_leg_thigh = Rototranslation_Rz(-length_leg,0,0, theta_knee);

R_0_leg = R_0_foot * R_foot_leg ; 
R_0_thigh = R_0_leg * R_leg_thigh ;

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

umt_length = sqrt(sum((Insertion - Origin).^2));
UMTLength = Function('UMTLength', ...
    {q, known_parameters}, {umt_length}, ...
    {'q', 'known_parameters'}, {'umt_length'});

momentArm = jacobian(umt_length, q);
getMomentArm = Function('MomentArm', ...
    {q, known_parameters}, {momentArm}, ...
    {'q', 'known_parameters'}, {'moment_arm'});



%% Activation Dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Groote F De, Kinney AL, Rao A V, Fregly BJ. Evaluation of Direct Collocation
% Optimal Control Problem Formulations for Solving the Muscle Redundancy Problem.
% Ann Biomed Eng. 2016; 44: 2922–2936. https://doi.org/10.1007/s10439-016-1591-9
% PMID: 27001399

% muscle activation is described by two nonlinear, first order differential
% equetaion

% parameter value
% taua = 0.015
% taud = 0.060
% b = 0.1

taua = SX.sym('Activation_time_constant',3);
taud = SX.sym('Desactivation_time_constant',3);
b = SX.sym('Transition_smoothness',3);

fa = 0.5.*tanh(b.*(e-a)) ;
da_dt = ((1 ./ taua .* (0.5 + 1.5 .* a)) .* (fa + 0.5) + ...
    ((0.5 + 1.5 .* a) ./ (taud)) .* (-fa + 0.5)) .* ...
    (e-a) ;


%% Muscle Contraction Dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Groote F De, Kinney AL, Rao A V, Fregly BJ. Evaluation of Direct Collocation
% Optimal Control Problem Formulations for Solving the Muscle Redundancy Problem.
% Ann Biomed Eng. 2016; 44: 2922–2936. https://doi.org/10.1007/s10439-016-1591-9
% PMID: 27001399

% muscle tendon properties to be identified
l0m = SX.sym('Optimal_fiber_length',3);
phi0 = SX.sym('Pennation_angle_at_muscle_optimal_fiber_length',3);
f0m = SX.sym('Maximal_isometric_muscle_force',3);
lst = SX.sym('Tendon_slack_length',3);


muscleTendonParameters =  vertcat(l0m,phi0,f0m,lst) ;


%% architectural parameter
normalizedMuscleLength = SX.sym('Normalized_muscle_length',3) ;
tendonLength = SX.sym('tendon_length',3) ;
muscleLength = SX.sym('muscle_length',3) ;


normalizedTendonLength = tendonLength ./ lst ; 



%% Tendon
% Tendon force-length (S1)
kT = 35 ; c1 = 0.200 ; c2 = 0.995 ; c3 = 0.250 ; % tendon parameters
normalizedTendonForce = c1 .* exp(kT .* (normalizedTendonLength - c2)) - c3 ; % tendon force function


tendonForce= normalizedTendonForce .* f0m ; 

getTendonForce = Function('getTendonForce', ...
    {q,known_parameters, muscleTendonParameters, normalizedTendonLength}, {tendonForce}, ...
    {'q','known_parameters','muscle_tendon_parameters','length_tendon_normalized'}, {'TendonForce'}) ;

%% Muscle
% Active muscle force-length (S2)
b11 = 0.815 ; b21 = 1.055 ; b31 = 0.162 ; b41 = 0.063 ; % first Gaussian coefficents
b12 = 0.433 ; b22 = 0.717 ; b32 = -0.030 ; b42 = 0.200 ; % second Gaussian coefficents
b13 = 0.100 ; b23 = 1.000 ; b33 = 0.354 ;  b43 = 0.000 ; % third Gaussian coefficents


NormalizedMuscleActiveForceLength = (b11 .* exp((-0.5.* (normalizedMuscleLength - b21).^2)./ (b31 + b41 .* normalizedMuscleLength))) + ...
    (b12 .* exp((-0.5.* (normalizedMuscleLength - b22).^2)./ (b32 + b42 .* normalizedMuscleLength))) + ...
    (b13 .* exp((-0.5.* (normalizedMuscleLength - b23).^2)./ (b33 + b43 .* normalizedMuscleLength))) ;



Normalized_Muscle_Active_Force_Length = Function('Normalized_Muscle_Active_Force_Length', {q,known_parameters,normalizedMuscleLength}, {NormalizedMuscleActiveForceLength}, ...
    {'q','known_parameters','length_muscle_normalized'}, {'NormalizedMuscleActiveForceLength'}) ;


MuscleActiveForceLength = NormalizedMuscleActiveForceLength .* f0m ; 

Muscle_Active_Force_Length = Function('Muscle_Active_Force_Length', ...
    {q,known_parameters,muscleTendonParameters,normalizedMuscleLength,a}, {MuscleActiveForceLength}, ...
    {'q','known_parameters','muscle_tendon_parameters','length_muscle_normalized','a'}, {'MuscleActiveForceLength'}) ;

% Passive muscle force-length (S3)
kpe = 4.0 ; e0 = 0.6 ;
NormalizedMusclePassiveForceLength = (exp(((kpe .* (normalizedMuscleLength - 1))./e0)) - 1)./ (exp(kpe) - 1) ;

Normalized_Muscle_Passive_Force_Length = Function('Normalized_Muscle_Passive_Force_Length', ...
    {q,known_parameters,normalizedMuscleLength}, {NormalizedMusclePassiveForceLength}, ...
    {'q','known_parameters','length_muscle_normalized'}, {'NormalizedMusclePassiveForceLength'}) ;

MusclePassiveForceLength = NormalizedMusclePassiveForceLength .* f0m ; 
Muscle_Passive_Force_Length = Function('Muscle_Passive_Force_Length', ...
    {q,known_parameters, muscleTendonParameters,normalizedMuscleLength}, {MusclePassiveForceLength}, ...
    {'q','known_parameters','muscle_tendon_parameters','length_muscle_normalized'}, {'MusclePassiveForceLength'}) ;

% Muscle force-velocity (S4)
% d1 -0.318
% d2 -8.149
% d3 -0.374
% d4 0.886
NormalizedMuscleForceVelocity = 1 ; % vitesse = 0


%% Forces function

normalizedMuscleForce = a .* NormalizedMuscleActiveForceLength .* NormalizedMuscleForceVelocity + NormalizedMusclePassiveForceLength ;
   
getNormalizedMuscleForce = Function('getNormalizedMuscleForce', ...
    {q, known_parameters, normalizedMuscleLength, a}, {normalizedMuscleForce}, ...
    {'q','known_parameters','length_muscle_normalized','a'}, {'NormalizedMuscleForce'}) ;

muscleForce = normalizedMuscleForce .* f0m ; 

getMuscleForce = Function('getMuscleForce', ...
    {q, known_parameters, muscleTendonParameters, normalizedMuscleLength, a}, {muscleForce}, ...
    {'q','known_parameters','muscle_tendon_parameters','length_muscle_normalized','a'}, {'MuscleForce'}) ;





%% Create all muscle-tendon functions
normalizeTendonLength = Function('normalizeTendonLength', ...
    {tendonLength, muscleTendonParameters}, {normalizedTendonLength}, ...
    {'length_tendon_normalized','muscle_tendon_parameters'}, {'length_tendon'});

muscleLengthNormalized = muscleLength ./ l0m ; 
normalizeMuscleLength = Function('normalizeMuscleLength', ...
    {muscleLength, muscleTendonParameters}, {muscleLengthNormalized}, ...
    {'muscleLength','muscleTendonParameters'}, {'normalizeMuscleLength'});

pennation_angle =  asin(l0m .*  sin(phi0) ./ normalizedMuscleLength) ;
Pennation_Angle = Function('Pennation_Angle', ...
    {normalizedMuscleLength, muscleTendonParameters}, {pennation_angle}, ...
    {'length_muscle_normalized','muscle_tendon_parameters'}, {'pennation_angle'});

normalizeTendonForce = Function('normalizeTendonForce', ...
    {q,known_parameters, normalizedTendonLength}, {normalizedTendonForce}, ...
    {'q','known_parameters','length_tendon_normalized'}, {'NormalizedTendonForce'}) ;









%% Computing Joint Moments and Angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Moment articulaire Mj(θ, t) = ∑ i=1 m (ri(θ) ⋅ Fimt(θ, t))
joint_moment = momentArm' * tendonForce ;
f_joint_moment = Function('joint_moment', {q,known_parameters,muscleTendonParameters,normalizedTendonLength}, {joint_moment}, ...
    {'q','known_parameters','muscle_tendon_parameters','length_tendon_normalized'}, {'joint_moment'}) ;


%% Muscle-tendon equilibrium
% determine muscle length such that TendonForce - (cos(pennation_angle) .* MuscleForce) = 0

length_tendon_normalized_num = normalizeTendonLength(tendonLength, muscle_tendon_parameters_num) ; 
length_muscle_normalized_num = normalizeMuscleLength(muscleLength, muscle_tendon_parameters_num) ; 

Ttest = getTendonForce(q_num,known_parameters_num,muscle_tendon_parameters_num,length_tendon_normalized_num) ;


g0 = tendonForce - (cos(pennation_angle) .* muscleForce) ; 
g1 = umt_length' - (cos(pennation_angle) .* muscleLength + tendonLength) ; 
g = Function('g',[tendonLength,muscleLength],[g0,g1]) ; 

equilibrateMuscleTendon = rootfinder('equilibrateMuscleTendon','newton',g) ; 




%% function d'optimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

unknown_parameters = horzcat(l0m,phi0,f0m,lst) ;
%%
% f minimisation quadratique
T = SX.sym('Mesured_torque',1) ;

f = (T -sum(joint_moment)).^2 ;
% g = Ft - Fm*cos(pennation_angle)

% %%
% nlp = struct('x',unknown_parameters, 'f',f, 'g',g);
% S = nlpsol('S', 'ipopt', nlp);
% sol = S('x0',x0, 'lbx',min_param,'ubx',max_param, 'lbg',0,'ubg',0); %il faut que lbg soit un vecteur de 0 de la taille de la contrainte
% Param = sol.x;
% disp(x_opt)