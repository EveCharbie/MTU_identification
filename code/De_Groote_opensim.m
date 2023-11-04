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


%% Topology 
% B : Bones
% J : Joins 
% df : degree of freedom
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total df = 6 df
% --> J : Hip (4 df : Rotx= 0 ; Roty = 0 ; Rotz = theta_hip ; X = x ; Y = y; Z = z )
%   --> B : Thigh
%    --> J : knee (1 df : Rotx= 0 ; Roty = 0 ; Rotz = theta_knee ; X = 0 ; Y = 0 ; Z = 0 )
%       --> B : Leg 
%           --> J : Ankle (1 dF: Rotx= 0 ; Roty = 0 ; Rotz = theta_ankle ; X = 0 ; Y = 0 ; Z = 0 )
%                --> B : Subtalar 
%                   --> J : talo_calc : talo_calcanéenne (0 dF: Rotx= 0 ; Roty = 0 ; Rotz = 0 ; X = 0 ; Y = 0 ; Z = 0 )
%                        --> B : Calc 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Musculoskeletal Geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% States
x = SX.sym('x');
y = SX.sym('y');
z = SX.sym('z');
theta_hip = SX.sym('theta_hip');
theta_knee = SX.sym('theta_knee');
theta_ankle = SX.sym('theta_ankle');
trans_talo_calc = SX.sym('translation_talo_calc',3);


q = vertcat(x,y,z,theta_hip, theta_knee, theta_ankle );

nMuscles = 3 ;
nBones = 4 ; 

%% Model parameters that are personalized by scaling
thigh = SX.sym('thigh',3);
leg = SX.sym('leg',3);
talus = SX.sym('talus',3);
clac = SX.sym('foot',3);

segment_geometry = vertcat(thigh, leg, talus ,clac);


%% rototanslation matrix
R_0_thigh = Rototranslation_Rz(x,y,z, theta_hip);
R_thigh_leg = Rototranslation_Rz(thigh(1),thigh(2),thigh(3), -theta_knee);
R_leg_talus = Rototranslation_Rz(leg(1),leg(2),leg(3), -theta_ankle );
R_talus_calc = Rototranslation_Rz(talus(1),talus(2),talus(3), 0 );


R_0_leg = R_0_thigh * R_thigh_leg ; 
R_0_talus = R_0_leg * R_leg_talus ; 
R_0_calc = R_0_talus * R_talus_calc ; 

%% Joint center
HJC = R_0_thigh(1:3, 4); % Hip
KJC = R_0_leg(1:3, 4); % Knee
AJC = R_0_talus(1:3, 4); % Knee
TJC = Rototranslate(R_0_calc, clac) ; % toe
CALC = R_0_calc(1:3, 4);  %calcaneum

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

known_parameters = vertcat(segment_geometry, muscle_insertion);

%% muscle origins and insertion in R0
Origin_tibialis_anterior = Rototranslate(R_0_leg, Local_Origin_tibialis_anterior);  
Insertion_tibialis_anterior = Rototranslate(R_0_calc, Local_Insertion_tibialis_anterior);
Origin_soleus = Rototranslate(R_0_leg,  Local_Origin_soleus) ;
Insertion_soleus = Rototranslate(R_0_calc, Local_Insertion_soleus);
Origin_gastrocnemius = Rototranslate(R_0_thigh, Local_Origin_gastrocnemius) ;
Insertion_gastrocnemius = Rototranslate(R_0_calc, Local_Insertion_gastrocnemius) ;

Origin = horzcat(Origin_tibialis_anterior, Origin_soleus, Origin_gastrocnemius);
Insertion = horzcat(Insertion_tibialis_anterior, Insertion_soleus, Insertion_gastrocnemius);
Markers = horzcat(HJC, KJC, AJC, TJC,CALC); 
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


