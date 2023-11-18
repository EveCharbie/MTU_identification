    %% test of the root finder [value NOT OK]
%Model_Opensim ; 
De_Groote_opensim ;

[known_parameters_num,muscle_tendon_parameters_num] = Opensim_extraction() ;

% input :  muscle activation and joint orientation  (random)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_num = [1, 1, 1]; 
q_num = [0, 0 , 0, 0, (0/180)*pi , (30/180)*pi ] ; 

% geometry 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[muscle_origin_in0, muscle_insertion_in0,muscle_viaPoint_in0, joint_centers] = ForwardKinematics(q_num, known_parameters_num);
umt_length = getUMTLength(q_num, known_parameters_num) ;
moment_arms = getMomentArm(q_num, known_parameters_num) ; 

%plotmodel(muscle_origin_in0, muscle_insertion_in0, joint_centers)
    



% muscle tendon parameter
%muscle_tendon_parameters_num =  [l0m_num , phi0_num , f0m_num , lst_num ] ; 
    % musculoskeletal states (q, bones length , muscletendon insertion) 
musculoskeletal_states_num = [q_num , known_parameters_num] ; 

p_num = horzcat(a_num,musculoskeletal_states_num , muscle_tendon_parameters_num) ;


x0 = [0.001, 0.001, 0.001 ,muscle_tendon_parameters_num(1:3)] ;  % tendonLengthening, fibre length   (TA SOL GAST)
x_num = equilibrateMuscleTendon(x0, p_num) ;

tendonLengthening_num = x_num(1:nMuscles)';
fiberLength_num = x_num(nMuscles+1:end)';


all_states_num = [musculoskeletal_states_num,  fiberLength_num, tendonLengthening_num , a_num] ;


Force_length_representation(all_states_num,muscle_tendon_parameters_num,casadiFun,NameMuscles)

