%% hypothetical data
% import model
De_Groote_opensim ;

[known_parameters_num,muscle_tendon_parameters_num] = Opensim_extraction() ;


%% trials
results = [] ;
header = {'Torque','q1','q2',...
    'activation_tibialis','activation_soleus','activation_gastrocnemius',...
    'fiber_tibialis','fiber_soleus','fiber_gastrocnemius',...
    'phi_tibialis','phi_soleus','phi_gastrocnemius'} ;

% input : Musculo skeletical configuration during trial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qankle = linspace(-40,10, 10)/180*pi ;
qknee = [0, 20, 40, 60, 80, 100]/180*pi ;

% input :  muscle activation  (random)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_num = [0, 0, 0; ...
    0, 1, 1;...
    1, 0, 0;];

ntrials = 0 ;
Data = [];

for i = 1 : size(a_num,1) % activation muscle
    for ii = 1 : length(qknee)
        for iii = 1 :  length(qankle)
            q_num = [0, 0 , 0, 0, qknee(ii), qankle(iii) ] ; %
            musculoskeletal_states_num = [q_num , known_parameters_num] ;
            neuromusculoskeletal_states_num = [a_num(i,:), musculoskeletal_states_num] ; 

            p_num = horzcat( neuromusculoskeletal_states_num, muscle_tendon_parameters_num) ;

            x0 = [0.001, 0.001, 0.001 ,muscle_tendon_parameters_num(1:3)] ;  % tendonLengthening, fibre length   (TA SOL GAST)
            x_num = full(equilibrateMuscleTendon(x0, p_num)) ; % x find


            tendonLengthening_num = x_num(1:nMuscles)';
            fiberLength_num = x_num(nMuscles+1:end)';
            rootedvariables_num = [fiberLength_num ,tendonLengthening_num] ; 
            all_states_num = [neuromusculoskeletal_states_num, rootedvariables_num] ;


            temp = full(casadiFun.getJointMoment(all_states_num,muscle_tendon_parameters_num)) ;
            Torque = temp(end) ;                                               % torque
            fLength_num = full(fiberLength_num) ;
            phi_num = full(casadiFun.getPennationAngle(all_states_num,muscle_tendon_parameters_num))' ;

            % output
            ntrials = ntrials +1 ;
            Data(ntrials,:) = [Torque, qknee(ii) , qankle(iii),  a_num(i,:),...
                fLength_num, phi_num, tendonLengthening_num ] ; % variables mesured 
            nbnan(ntrials) = sum(sum(isnan(Data(ntrials,:)))) ;
        end
    end
end


fprintf('Number of nan for the trial :\n')
fprintf(num2str(nbnan))
fprintf('\n')