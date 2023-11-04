    %% test of the root finder 
Model_Opensim ; 

known_parameters_num = [segment_geometry_num; muscle_origin; muscle_insersion];
q_num = [0, 0 , 0, 0, (10/180)*pi , (0/180)*pi ] ; 
% geometry 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[muscle_origin_in0, muscle_insertion_in0, joint_centers] = ForwardKinematics(q_num, known_parameters_num);
umt_length = getUMTLength(q_num, known_parameters_num) ; 
moment_arms = getMomentArm(q_num, known_parameters_num) ; 

plotmodel(muscle_origin_in0, muscle_insertion_in0, joint_centers)
    
% neuro musculo skeltical states (random)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_num = [.2, .3, .4]; 
  

% muscle tendon parameter
muscle_tendon_parameters_num =  [l0m_num , phi0_num , f0m_num , lst_num ] ; 
    % musculoskeletal states (q, bones length , muscletendon insertion) 
musculoskeletal_states_num = [q_num , known_parameters_num'] ; 

p_num = horzcat(a_num,musculoskeletal_states_num , muscle_tendon_parameters_num) ;


x0 = [0.001, 0.001, 0.001 ,l0m_num] ;  % tendonLengthening, fibre length   (TA SOL GAST)
equilibrateMuscleTendon(x0, p_num)