%% hypothetical data
% import model
De_Groote_opensim ;

[known_parameters_num,muscle_tendon_parameters_num] = Opensim_extraction() ;


%% trials
nTrials = 6 ;

results = {} ;
headler = {'Torque','q1','q2',...
    'activation_tibialis','activation_soleus','activation_gastrocnemius',...
    'momentarm_tibialis','momentarm_soleus','momentarm_gastrocnemius',...
    'UMT_tibialis','UMT_soleus','UMT_gastrocnemius',...
    'fiber_tibialis','fiber_soleus','fiber_gastrocnemius',...
    'phi_tibialis','phi_soleus','phi_gastrocnemius',...
     'tendon_tibialis','tendon_soleus','tendon_gastrocnemius'} ;

% Musculo skeletical configuration during trial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qankle = linspace(-50,20, 100) ; 
qknee = [0, 20, 40, 60, 80, 100] ;

% input :  muscle activation and joint orientation  (random)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_num = [1, 1, 1];

for i = 1 : nTrials
    for ii = 1 :  length(qankle)
        q_num = [0, 0 , 0, 0, (qknee(i)/180)*pi , (qankle(ii)/180)*pi ] ; %
        musculoskeletal_states_num = [q_num , known_parameters_num] ;

        p_num = horzcat(a_num,musculoskeletal_states_num , muscle_tendon_parameters_num) ;

        x0 = [0.001, 0.001, 0.001 ,muscle_tendon_parameters_num(1:3)] ;  % tendonLengthening, fibre length   (TA SOL GAST)
        x_num = equilibrateMuscleTendon(x0, p_num) ; % x find
        

        %if sum(isnan(full(x_num)))== 0
        tendonLengthening_num = x_num(1:nMuscles)';
        fiberLength_num = x_num(nMuscles+1:end)';
        all_states_num = [musculoskeletal_states_num,  fiberLength_num, tendonLengthening_num , a_num] ;
        

        temp = full(casadiFun.getJointMoment(all_states_num,muscle_tendon_parameters_num)) ; 
        Torque = temp(end) ;                                           % torque
        %Q1
        %Q2
        %Activation
        temp = full(casadiFun.getMomentArm(q_num,known_parameters_num));
        ma = temp(:,end)' ;                                                % moment arm
        UMTlength_num = full(casadiFun.getUMTLength(q_num,known_parameters_num))' ; % UMT length
        fLength_num = full(fiberLength_num) ;
        phi_num = full(casadiFun.getPennationAngle(all_states_num,muscle_tendon_parameters_num))' ;
        tendonlength_num = full(casadiFun.normalizeTendonLength(all_states_num,muscle_tendon_parameters_num))' .* muscle_tendon_parameters_num(10:12) ; 

        Data(ii,:) = [Torque, qknee(i) , qankle(ii), a_num, ma,...
            UMTlength_num, fLength_num, ...
            phi_num, tendonlength_num] ; 

    end
        results{i} = Data ; 
        nbnan(i) = sum(sum(isnan(Data))) ;
        Data = [] ;
end
fprintf('Number of nan for the trial :\n')
fprintf(num2str(nbnan))
fprintf('\n')