%% test rootfinder of muscle tendon equilibrium

example_state_param; 


% x = vertcat(tendonLengthening, fiberLength);
% param_ = vertcat(states, known_parameters, muscleTendonParameters);

param_ = [state_num; known_parameters_num; muscle_tendon_parameters_num];

x_ = equilibrateMuscleTendon([0, 0, 0, .1, .1, .1], param_);  
tendonLengthening_ = x_(1:nMuscles);
fiberLength_ = x_(nMuscles+1:end);

all_states_num = [state_num; tendonLengthening_; fiberLength_];

tendonForce_ = getTendonForce(all_states_num, known_parameters_num, muscle_tendon_parameters_num); 
muscleForce_ = getMuscleForce(all_states_num, known_parameters_num, muscle_tendon_parameters_num); 
muscleAForce_ = getMuscleActiveForce(all_states_num, known_parameters_num, muscle_tendon_parameters_num);
musclePForce_ = getMusclePassiveForce(all_states_num, known_parameters_num, muscle_tendon_parameters_num);



% 
% getPennationAngle = Function('getPennationAngle', ...
%     {all_states, muscleTendonParameters}, {pennationAngle}, ...
%     {'all_states','muscle_tendon_parameters'}, {'pennation_angle'});
% 
% normalizeTendonForce = Function('normalizeTendonForce', ...
%     {all_states, known_parameters, muscleTendonParameters}, {normalizedTendonForce}, ...
%     {'all_states','known_parameters', 'muscleTendonParameters'}, {'NormalizedTendonForce'});
% 
% getMusclePassiveForce = Function('getMusclePassiveForce', ...
%     {all_states, known_parameters, muscleTendonParameters}, {musclePassiveForce}, ...
%     {'all_states','known_parameters','muscle_tendon_parameters'}, {'MusclePassiveForce'}) ;
% 
% getMuscleActiveForce = Function('getMuscleActiveForce', ...
%     {all_states, known_parameters, muscleTendonParameters}, {MuscleActiveForceLength}, ...
%     {'all_states','known_parameters','muscle_tendon_parameters'}, {'MuscleActiveForce'}) ;
% 
% 
