% test casadi functions
%run('.\De_Groote'); 
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


[muscle_origin_in0, muscle_insertion_in0, joint_centers] = ForwardKinematics(q_num, known_parameters_num);
umt_length = UMTLength(q_num, known_parameters_num);
moment_arms = MomentArm(q_num, known_parameters_num);


plotmodel(muscle_origin_in0, muscle_insertion_in0, joint_centers)
