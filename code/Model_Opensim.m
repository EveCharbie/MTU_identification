% New Model Opensim 

% import the model (De_Groot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
De_Groote_opensim ; 

% Muscle insetion extract from  Opensim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X : 
% Y : longitudial axis 
% foot; leg; thigh

% position knee dans leg  
% kne location_in_parent (thigh)
knee = [0, -0.39, 0] ; 

% ankle location_in_parent (leg) 
ankle = [0, -0.43, 0] ; 

% subtalar location_in_parent (talus) 
subtalar = [-0.04877, -0.04195, 0 ]; %  0.00792

% mtp location_in_parent (calc) 
toe = [0.1788, -0.002, 0]; % 0.00108

segment_geometry_num = [knee'; ankle'; subtalar'; toe'];

%    local TA SOL GAST (X, Y Z)
    % Origin 
TA_origin = [0.018, -0.162, 0] ; % z = -0.011 non 0
SOL_origin = [-0.002, -0.153, 0] ;  % z = -0.007
GASTm_origin = [-0.019, -0.393, 0] ; % z =0.024
GASTl_origin = [-0.022, -0.395, 0] ; % z = -0.027
GAST_origin = mean([GASTm_origin; GASTl_origin]) ; % z = -0.027

    % via points
TA_viaPoint = [0.033, -0.395, 0] ; % z = 0.018 non 0 --> referentiel Tibia
GASTm_ViaPoint = [-0.03, -0.402, 0] ; % z= 0.026
GASTl_ViaPoint = [-0.03, -0.402, 0] ; % z= 0.026
GAST_ViaPoint =  mean([GASTm_ViaPoint; GASTl_ViaPoint]) ;


    % Insertions  
TA_insersion = [0.117, 0.018, 0] ; % z = 0.03 non 0
SOL_insersion = [0, 0.0310, 0] ;  % z =0.005
GASTm_insersion = [0, 0.031, 0] ; %z = 0.005
GASTl_insersion = [0, 0.031, 0] ; %z = 0.005
GAST_insersion = mean([GASTm_insersion; GASTl_insersion]) ;


%%
muscle_origin = [TA_origin' ; 
    SOL_origin' ; 
    GAST_origin' ];

muscle_insersion = [TA_insersion' ; 
    SOL_insersion' ; 
    GAST_insersion' ];




% % Muscle parameter extract from  Opensim
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TA
TA_l0m = 0.098 ; 
TA_phi0 = 0.08726646 ; 
TA_f0m = 905.0 ; 
TA_lst = 0.223 ; 

    % SOL
SOL_l0m = 0.05 ; 
SOL_phi0 = 0.43633231 ; 
SOL_f0m = 3549.0 ; 
SOL_lst = 0.25 ; 

    %GM (medial, lateral)
GAST_l0m = mean([0.06, 0.064]) ; 
GAST_phi0 = mean([0.29670597, 0.13962634]) ; 
GAST_f0m = sum([1558.0 ,683.0]); 
GAST_lst = mean([0.39, 0.38]); 

    % define muscle tendon parameter 
l0m_num = [TA_l0m ,  SOL_l0m , GAST_l0m ] ; 
phi0_num = [TA_phi0  , SOL_phi0 , GAST_phi0  ] ;
f0m_num = [TA_f0m , SOL_f0m , GAST_f0m] ; 
lst_num = [TA_lst , SOL_lst , GAST_lst ] ; 
