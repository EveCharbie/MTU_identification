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
Itad = SX.sym('Local_Insertion_tibialis_anterior',3);
Itap = SX.sym('Local_Origin_tibialis_anterior',3);
Isod = SX.sym('Local_Insertion_soleus',3);
Isop = SX.sym('Local_Origin_soleus',3);
Igad = SX.sym('Local_Insertion_gastrocnemius',3);
Igap = SX.sym('Local_Origin_gastrocnemius',3);
muscle_insersion = vertcat(Itad, Itap, Isod, Isop, Igad, Igap);

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
% lm = SX.sym('mesured_muscle_length',3);
% phi = SX.sym('muscle_pennation_angle',3);
a_t= SX.sym('normalized_neuronal activation',3);


%% model geometric 
% Rototranslation
% axe longitudinal x
R_0_foot = SX.eye(4);
R_0_foot(1:2, 1:2) = [ cos(theta_foot), - sin(theta_foot) ; ...
    sin(theta_foot), cos(theta_foot)] ;
R_0_foot(1:3,4) = [x, y, z]';

R_foot_leg = SX.eye(4);
R_foot_leg(1:2, 1:2) =[ cos(theta_ankle), - sin(theta_ankle); ...
    sin(theta_ankle), cos(theta_ankle)] ;
R_foot_leg(1,4) = length_foot;

R_leg_limb = SX.eye(4);
R_leg_limb(1:2, 1:2) = [ cos(theta_knee), - sin(theta_knee); ...
    sin(theta_knee), cos(theta_knee)] ;
R_leg_limb(1,4) = length_leg ;

R_0_leg = R_0_foot * R_foot_leg ; %T0_2 = T0_1 * T1_2
R_0_limb = R_0_leg * R_leg_limb ;

    %% Joint center 
% TJC = R_0_foot * [0;0;0;1] ; %toe 
% AJC = R_0_foot * [0;length_foot;0; 1] ; % ankle
% KJC = R_0_ankle * [0;length_leg;0; 1] ; % knee
% HJC = R_0_knee * [0;length_thigh;0; 1] ; % hip

TJC = R_0_foot(:,4) ; 
AJC = R_0_leg(:,4) ; % ankle
KJC = R_0_limb(:,4); % knee
HJC = R_0_limb * [length_thigh; 0; 0; 1] ; % hip

    %% muscle insertion et origine
%tibialis
ITAD = R_0_foot * [Itad; 1] ; 
ITAP = R_0_leg * [Itap; 1] ; 
%soleus
ISOD = R_0_foot * [Isod; 1] ; 
ISOP = R_0_leg * [Isop; 1] ; 
% gastrocnemius
IGAD = R_0_foot * [Igad; 1] ; 
IGAP = R_0_limb * [Igap; 1] ; 


Origin = horzcat(ITAP(1 :3),ISOP(1 :3),IGAP(1 :3));
Insertion = horzcat(ITAD(1 :3),ISOD(1 :3),IGAD(1 :3));
Markers = horzcat(TJC(1 :3), AJC(1 :3), KJC(1 :3), HJC(1 :3));


                        %% functions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ForwardKinematics = Function('ForwardKinematics', {q, known_parameters}, {Origin, Insertion, Markers}, ...
    {'q', 'known_parameters'}, {'Origin', 'Insertion', 'Markers'}) ; 

umt_length = sqrt(sum((Insertion - Origin).^2));
UMTLength = Function('UMTLength', {q, known_parameters}, {umt_length}, ...
    {'q', 'known_parameters'}, {'muscle_length'}); 

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
muscle_tendon_parameters =  horzcat(l0m,phi0,f0m,lst) ; 

lt = SX.sym('Normalized_tendon_length',3) ;
lm = SX.sym('Normalized_muscle_length',3) ;
phi = SX.sym('Pennation_angle',3) ;

% input a

%% Tendon 
% Tendon force-length (S1)
kT = 35 ; c1 = 0.200 ; c2 = 0.995 ; c3 = 0.250 ; % tendon parameters 
ft = c1 .* exp(kT .* (lt - c2)) - c3 ; % tendon force function 

%% Muscle 
% Active muscle force-length (S2)
b11 = 0.815 ; b21 = 1.055 ; b31 = 0.162 ; b41 = 0.063 ; % first Gaussian parameter
b12 = 0.433 ; b22 = 0.717 ; b32 = -0.030 ; b42 = 0.200 ; % second Gaussian parameter
b13 = 0.100 ; b23 = 1.000 ; b33 = 0.354 ;  b43 = 0.000 ; % third Gaussian parameter

fma = (b11 .* exp((-0.5.* (lm - b21).^2)./ (b31 + b41 .* lm))) + ...
    (b12 .* exp((-0.5.* (lm - b22).^2)./ (b32 + b42 .* lm))) + ...
    (b13 .* exp((-0.5.* (lm - b23).^2)./ (b33 + b43 .* lm))) ;

% Passive muscle force-length (S3)
        kpe = 4.0 ; e0 = 0.6 ; 
fmp = (exp(((kpe .* (lm - 1))./e0)) - 1)./ (exp(kpe) - 1) ; 

% Muscle force-velocity (S4)
        % d1 -0.318
        % d2 -8.149
        % d3 -0.374
        % d4 0.886
fva = 1 ; % vitesse = 0 


%% Forces function 
%tendon forces 
Ft = ft .* f0m ; 
FT = Function('FT', {q,known_parameters,muscle_tendon_parameters,lt}, {Ft}, ...
    {'q','known_parameters','muscle_tendon_parameters','lt'}, {'Ft'}) ; 

Fm = (a .* fma .* fva + fmp) .* f0m ; 
FM = Function('FM', {q,known_parameters,muscle_tendon_parameters,lm,a}, {Fm}, ...
    {'q','known_parameters','muscle_tendon_parameters','lm','a'}, {'Fm'}) ; 

                        %% Computing Joint Moments and Angles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Moment articulaire Mj(θ, t) = ∑ i=1 m (ri(θ) ⋅ Fimt(θ, t))
moment_articulaire = moment_arm .* Ft ; 
Momentarticualire = Function('Momentarticualire', {q,known_parameters,muscle_tendon_parameters,lt}, {moment_articulaire}, ...
    {'q','known_parameters','muscle_tendon_parameters','lt'}, {'moment_articulaire'}) ; 





                        %% function d'optimisation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 1. « Trouver »  des valeurs de paramètres qui font un peu de sens  :
% issues d'un scale opensim
% soleus 
Sol_f0m = 3600 ; 
Sol_l0m = 0.031 ; 
Sol_phi0 = 0.20943951 ; %rad
Sol_lst = 0.31 ; 

% Gast 
Gast_f0m = 2500 ; 
Gast_l0m = 0.09 ; 
Gast_phi0 = 0.29670597 ; %rad
Gast_lst = 0.36 ; 

% Ta 
Ta_f0m = 3000 ; 
Ta_l0m = 0.098 ; 
Ta_phi0 = 0.29670597 ; %rad
Ta_lst = 0.08726646 ; 

% 2. Déterminer des bornes physiologiques pour les paramètres (en faire des vecteurs : min_param et max_param)
%min_param = 
%Ta_phi0 = 0
% Ta_l0m = 0
%max_param = 
% Ta_phi0 = pi/2
% Ta_l0m = 0.25
% Ta_lst = 


% 3. Générer tes essais (en faire une liste et appeler les fonctions pour avoir les couples de référence et les longueurs du tendon et du muscle)

% 4. Ecrire le pb d’optimisation
unknown_parameters = horzcat(l0m,phi0,f0m,lst) ; 
%%
% f minimisation quadratique 
T = SX.sym('Mesured_torque',1) ;

f = (T -sum(moment_articulaire)).^2 ; 
g = Ft - Fm*cos(phi)

%%
nlp = struct('x',unknown_parameters, 'f',f, 'g',g);
S = nlpsol('S', 'ipopt', nlp);
sol = S('x0',x0, 'lbx',min_param,'ubx',max_param, 'lbg',0,'ubg',0); %il faut que lbg soit un vecteur de 0 de la taille de la contrainte
Param = sol.x;
disp(x_opt)