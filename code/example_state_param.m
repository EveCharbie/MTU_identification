%% test du model géometrique (OK)
De_Groote ; 
% foot; leg; thigh
segment_length_num = [.2; .4; .45];

% Origin TA, SO, GA; insertion TA; SO; GA
muscle_origin = [-.3; .01 ; 0 ; 
    -.3 ; -.01 ; 0 ; 
    -.02; -.01 ; 0 ];

muscle_insersion = [-.17; 0.01; 0 ; 
    -.23; -.02; 0 ; 
    -.23; -.02; 0];

known_parameters_num = [segment_length_num;muscle_origin;muscle_insersion];

q_num = [0; 0; 0; -20/180*pi; 10/180*pi; 50/180*pi];
activation_ = [.2; .3; .4]; 
state_num = [q_num; activation_];

% TA; SO; GA;
l0m_num = [0.098 ;  0.031 ; 0.09] ; 
phi0_num = [0.29670597 ; 0.20943951 ; 0.29670597 ] ;
f0m_num = [3000 ; 3600 ; 2500] ; 
lst_num = [0.08726646 ; 0.31 ; 0.36 ] ; 
muscle_tendon_parameters_num =  [l0m_num ; phi0_num ; f0m_num ; lst_num ] ; 
