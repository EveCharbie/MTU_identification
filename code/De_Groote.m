%% Casadi
addpath('C:\Users\Stage\Desktop\Neuromusculoskeletal Modeling\Casadi')
import casadi.*



%% Musculoskeletal Geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% definition du model
% Joint Center :
% TJC : toe joint center
% AJC : ankle joint center
% KJC : knee joint center
% HJC : hip joint center

% degre de liberté
x = SX.sym('x');
y = SX.sym('y');
z = SX.sym('z');
theta_foot = SX.sym('theta_foot');
theta_ankle = SX.sym('theta_ankle');
theta_knee = SX.sym('theta_knee');
q = vertcat(x,y,z,theta_foot, theta_ankle, theta_knee);

% Segment (bones)
% foot
% leg
% thigh

length_foot = SX.sym('length_foot');
length_leg = SX.sym('length_leg');
length_thigh = SX.sym('length_thigh');
segment_length = vertcat(length_foot, length_leg, length_thigh);

% Insersion muscle
% Ita : insersion distal du tibilais anterieur (sur le pied)
% Itap : insersion proximal du tibilais anterieur (sur le tibia_t)
% Isod : insersion distal du soléaire (sur le pied, calc)
% Isop : insersion proximal du soléaire (sur le tibia)
% Igad : insersion distal du gastrocnemius (sur le pied, calc)
% Igap : insersion proximal du gastrocnemius (sur la cuisse)

% position dans le repère local x,y, z (=0)
Local_Insertion_tibialis_anterior = SX.sym('Local_Insertion_tibialis_anterior',3);
Local_Origin_tibialis_anterior = SX.sym('Local_Origin_tibialis_anterior',3);
Local_Insertion_soleus = SX.sym('Local_Insertion_soleus',3);
Local_Origin_soleus = SX.sym('Local_Origin_soleus',3);
Local_Insertion_gastrocnemius = SX.sym('Local_Insertion_gastrocnemius',3);
Local_Origin_gastrocnemius = SX.sym('Local_Origin_gastrocnemius',3);
muscle_insersion = vertcat( Local_Origin_tibialis_anterior, ...
    Local_Origin_soleus, ...
    Local_Origin_gastrocnemius, ...
    Local_Insertion_tibialis_anterior, ...
    Local_Insertion_soleus, ...
    Local_Insertion_gastrocnemius );

known_parameters = vertcat(segment_length, muscle_insersion) ;

% Estimated parameter
% Fom : maximal isometric muscle force (Fom),
% lom  : optimal fiber length (ℓom),
% lst : tendon slack length (ℓst),
% alpha0 : optimal fiber length(φo).
% estimated_parameters = vertcat(Fom, lom,lst,alpha0);

%     % parmètre musculaire
% lambda = SX.sym('Percentage_change_in_optimal_fibre_length',3);
%
% % mesure
% length_muscle_normalized = SX.sym('mesured_muscle_length',3);
% pennation_angle = SX.sym('muscle_pennation_angle',3);
a_t= SX.sym('normalized_neuronal activation',3);


%% model geometric
% Rototranslation
% axe longitudinal x

R_0_foot = Rototranslation_Rz(x,y,z, theta_foot);
R_foot_leg = Rototranslation_Rz(-length_foot,0,0, theta_ankle-pi/2);
R_leg_thigh = Rototranslation_Rz(-length_leg,0,0, theta_knee);

R_0_leg = R_0_foot * R_foot_leg ; 
R_0_thigh = R_0_leg * R_leg_thigh ;

%% Joint center
% TJC = R_0_foot * [0;0;0;1] ; %toe
% AJC = R_0_foot * [0;length_foot;0; 1] ; % ankle
% KJC = R_0_ankle * [0;length_leg;0; 1] ; % knee
% HJC = R_0_knee * [0;length_thigh;0; 1] ; % hip

TJC = R_0_foot(1:3, 4);
AJC = R_0_leg(1:3, 4); % ankle
KJC = R_0_thigh(1:3, 4); % knee
HJC = Rototranslate(R_0_thigh,  [-length_thigh; 0; 0;]) ; % hip

%% muscle insertion et origine
%tibialis
Origin_tibialis_anterior = Rototranslate(R_0_leg, Local_Origin_tibialis_anterior);  
Insertion_tibialis_anterior = Rototranslate(R_0_foot, Local_Insertion_tibialis_anterior);
%soleus
Origin_soleus = Rototranslate(R_0_leg,  Local_Origin_soleus) ;
Insertion_soleus = Rototranslate(R_0_foot, Local_Insertion_soleus);
% gastrocnemius
Origin_gastrocnemius = Rototranslate(R_0_thigh, Local_Origin_gastrocnemius) ;
Insertion_gastrocnemius = Rototranslate(R_0_foot, Local_Insertion_gastrocnemius) ;


Origin = horzcat(Origin_tibialis_anterior, Origin_soleus, Origin_gastrocnemius);
Insertion = horzcat(Insertion_tibialis_anterior, Insertion_soleus, Insertion_gastrocnemius);
Markers = horzcat(TJC, AJC, KJC, HJC);


%% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ForwardKinematics = Function('ForwardKinematics', {q, known_parameters}, {Origin, Insertion, Markers}, ...
    {'q', 'known_parameters'}, {'Origin', 'Insertion', 'Markers'}) ;

umt_length = sqrt(sum((Insertion - Origin).^2));
UMTLength = Function('UMTLength', {q, known_parameters}, {umt_length}, ...
    {'q', 'known_parameters'}, {'umt_length'});

moment_arm = jacobian(umt_length, q);
MomentArm = Function('MomentArm', {q, known_parameters}, {moment_arm}, ...
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

e = SX.sym('Neuromuscular_activation',3);
a = SX.sym('Muscle_activation',3);
taua = SX.sym('Activation_time_constante',3);
taud = SX.sym('Desactivation_time_constante',3);
b = SX.sym('Parameter_determining_transition_smoothness',3);

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

% muscle tendon properties
l0m = SX.sym('Optimal_fiber_length',3) ;
phi0 = SX.sym('Pennation_angle_at_muscle_optimal_fiber_length',3) ;
f0m = SX.sym('Maximal_isometric_muscle_force',3) ;
lst = SX.sym('Tendon_slack_length',3) ;
muscle_tendon_parameters =  vertcat(l0m,phi0,f0m,lst) ;


%% architectural parameter
length_tendon_normalized = SX.sym('Normalized_tendon_length',3) ;
length_muscle_normalized = SX.sym('Normalized_muscle_length',3) ;
length_tendon = SX.sym('tendon_length',3) ;
length_muscle = SX.sym('muscle_length',3) ;


LengthTendonNomralized = length_tendon ./ lst ; 
Length_Tendon_Nomralized = Function('Length_Tendon_Nomralized', {length_tendon,muscle_tendon_parameters}, {LengthTendonNomralized}, ...
    {'length_tendon_normalized','muscle_tendon_parameters'}, {'length_tendon'}) ;

LengthMuscleNomralized = length_muscle ./ l0m ; 
Length_Muscle_Nomralized = Function('Length_Muscle_Nomralized', {length_muscle,muscle_tendon_parameters}, {LengthMuscleNomralized}, ...
    {'length_muscle_normalized','muscle_tendon_parameters'}, {'LengthMuscleNomralized'}) ;

pennation_angle =  asin(l0m .*  sin(phi0) ./ length_muscle_normalized) ;
Pennation_Angle = Function('Pennation_Angle', {length_muscle_normalized,muscle_tendon_parameters}, {pennation_angle}, ...
    {'length_muscle_normalized','muscle_tendon_parameters'}, {'pennation_angle'}) ;

%% Tendon
% Tendon force-length (S1)
kT = 35 ; c1 = 0.200 ; c2 = 0.995 ; c3 = 0.250 ; % tendon parameters
NormalizedTendonForce = c1 .* exp(kT .* (length_tendon_normalized - c2)) - c3 ; % tendon force function

Normalized_Tendon_Force = Function('Normalized_Tendon_Force', {q,known_parameters,length_tendon_normalized}, {NormalizedTendonForce}, ...
    {'q','known_parameters','length_tendon_normalized'}, {'NormalizedTendonForce'}) ;

TendonForce= NormalizedTendonForce .* f0m ; 

Tendon_Force = Function('Tendon_Force', {q,known_parameters,muscle_tendon_parameters,length_tendon_normalized}, {TendonForce}, ...
    {'q','known_parameters','muscle_tendon_parameters','length_tendon_normalized'}, {'TendonForce'}) ;

%% Muscle
% Active muscle force-length (S2)
b11 = 0.815 ; b21 = 1.055 ; b31 = 0.162 ; b41 = 0.063 ; % first Gaussian coefficents
b12 = 0.433 ; b22 = 0.717 ; b32 = -0.030 ; b42 = 0.200 ; % second Gaussian coefficents
b13 = 0.100 ; b23 = 1.000 ; b33 = 0.354 ;  b43 = 0.000 ; % third Gaussian coefficents


NormalizedMuscleActiveForceLength = (b11 .* exp((-0.5.* (length_muscle_normalized - b21).^2)./ (b31 + b41 .* length_muscle_normalized))) + ...
    (b12 .* exp((-0.5.* (length_muscle_normalized - b22).^2)./ (b32 + b42 .* length_muscle_normalized))) + ...
    (b13 .* exp((-0.5.* (length_muscle_normalized - b23).^2)./ (b33 + b43 .* length_muscle_normalized))) ;



Normalized_Muscle_Active_Force_Length = Function('Normalized_Muscle_Active_Force_Length', {q,known_parameters,length_muscle_normalized}, {NormalizedMuscleActiveForceLength}, ...
    {'q','known_parameters','length_muscle_normalized'}, {'NormalizedMuscleActiveForceLength'}) ;


MuscleActiveForceLength = NormalizedMuscleActiveForceLength .* f0m ; 

Muscle_Active_Force_Length = Function('Muscle_Active_Force_Length', {q,known_parameters,muscle_tendon_parameters,length_muscle_normalized,a}, {MuscleActiveForceLength}, ...
    {'q','known_parameters','muscle_tendon_parameters','length_muscle_normalized','a'}, {'MuscleActiveForceLength'}) ;

% Passive muscle force-length (S3)
kpe = 4.0 ; e0 = 0.6 ;
NormalizedMusclePassiveForceLength = (exp(((kpe .* (length_muscle_normalized - 1))./e0)) - 1)./ (exp(kpe) - 1) ;

Normalized_Muscle_Passive_Force_Length = Function('Normalized_Muscle_Passive_Force_Length', {q,known_parameters,length_muscle_normalized}, {NormalizedMusclePassiveForceLength}, ...
    {'q','known_parameters','length_muscle_normalized'}, {'NormalizedMusclePassiveForceLength'}) ;

MusclePassiveForceLength = NormalizedMusclePassiveForceLength .* f0m ; 
Muscle_Passive_Force_Length = Function('Muscle_Passive_Force_Length', {q,known_parameters,muscle_tendon_parameters,length_muscle_normalized}, {MusclePassiveForceLength}, ...
    {'q','known_parameters','muscle_tendon_parameters','length_muscle_normalized'}, {'MusclePassiveForceLength'}) ;

% Muscle force-velocity (S4)
% d1 -0.318
% d2 -8.149
% d3 -0.374
% d4 0.886
NormalizedMuscleForceVelocity = 1 ; % vitesse = 0


%% Forces function

NormalizedMuscleForce = a .* NormalizedMuscleActiveForceLength .* NormalizedMuscleForceVelocity + NormalizedMusclePassiveForceLength ;
   
Normalized_Muscle_Force = Function('Normalized_Muscle_Force', {q,known_parameters,length_muscle_normalized,a}, {NormalizedMuscleForce}, ...
    {'q','known_parameters','length_muscle_normalized','a'}, {'NormalizedMuscleForce'}) ;

MuscleForce = NormalizedMuscleForce .* f0m ; 

Muscle_Force = Function('Muscle_Force', {q,known_parameters,muscle_tendon_parameters,length_muscle_normalized,a}, {MuscleForce}, ...
    {'q','known_parameters','muscle_tendon_parameters','length_muscle_normalized','a'}, {'MuscleForce'}) ;

%% Computing Joint Moments and Angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Moment articulaire Mj(θ, t) = ∑ i=1 m (ri(θ) ⋅ Fimt(θ, t))
moment_articulaire = moment_arm' * TendonForce ;
Momentarticualire = Function('Momentarticualire', {q,known_parameters,muscle_tendon_parameters,length_tendon_normalized}, {moment_articulaire}, ...
    {'q','known_parameters','muscle_tendon_parameters','length_tendon_normalized'}, {'moment_articulaire'}) ;


%% function d'optimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

unknown_parameters = horzcat(l0m,phi0,f0m,lst) ;
%%
% f minimisation quadratique
T = SX.sym('Mesured_torque',1) ;

f = (T -sum(moment_articulaire)).^2 ;
% g = Ft - Fm*cos(pennation_angle)

% %%
% nlp = struct('x',unknown_parameters, 'f',f, 'g',g);
% S = nlpsol('S', 'ipopt', nlp);
% sol = S('x0',x0, 'lbx',min_param,'ubx',max_param, 'lbg',0,'ubg',0); %il faut que lbg soit un vecteur de 0 de la taille de la contrainte
% Param = sol.x;
% disp(x_opt)