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
qankle = [-30, -20, -10, 0, 10, 20, 30, ...
    25, 15, 5, -5, -15, -25]...
    /180*pi ;
qknee = [0, 20, 40, 60, 80, 70, 50, 30, 10]/180*pi ;

% input :  muscle activation  (random)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_num = [0, 0, 0; ...
    0, 0.2, 0.2;...
    0, 0.4, 0.4;...
    0, 0.6, 0.6;...
    0, 0.8, 0.8;...
    0, 1, 1] ;

ntrials = 0 ;
ntrialsfail = 0 ; 
Data = [];

%%
    % muscle tendon parameter 
muscle_tendon_parameters_ta = muscle_tendon_parameters_num([1,4,7,10]) ; 
muscle_tendon_parameters_sol = muscle_tendon_parameters_num([2,5,8,11]) ; 
muscle_tendon_parameters_gast = muscle_tendon_parameters_num([3,6,9,12]) ; 
    % inital guess
unknown_ta = [ones(1,3).* .0001 ,muscle_tendon_parameters_ta(1:2)] ;
unknown_sol = [ones(1,3).* .0001 ,muscle_tendon_parameters_sol(1:2)] ;
unknown_gast = [ones(1,3).* .0001 ,muscle_tendon_parameters_gast(1:2)] ;


    compt = 0 ; 
for i = 1 : size(a_num,1) % activation muscle
    for ii = 1 : length(qknee)
        for iii = 1 :  length(qankle)
            compt = compt+1 ; 
            % Neuromusculoskeletal configuration
            q_num = [0, 0 , 0, 0, qknee(ii), qankle(iii) ] ; 
            musculoskeletal_states_num = [q_num , known_parameters_num] ;
            neuromusculoskeletal_states_num = [a_num(i,:), musculoskeletal_states_num] ; 
            p_num = horzcat( neuromusculoskeletal_states_num, muscle_tendon_parameters_num) ;
            
            % Get UMT length 
            temp = casadiFun.getUMTLength(musculoskeletal_states_num) ; 
            length_UMT_ta(1,compt) = full( temp(1)) ;
            length_UMT_sol(1,compt) = full( temp(2)) ;
            length_UMT_gast(1,compt) = full( temp(3)) ;
            
            % known parameter 
            known_ta = [a_num(i,1), length_UMT_ta(1,compt) , muscle_tendon_parameters_ta] ;
            known_sol = [a_num(i,2), length_UMT_sol(1,compt) , muscle_tendon_parameters_sol] ;
            known_gast = [a_num(i,3), length_UMT_gast(1,compt) , muscle_tendon_parameters_gast] ;

            % Equilibrium
            unknown_ta = full(casadiFun.equilibrateMuscleTendonSingleMuscle(unknown_ta, known_ta)) ; % x find            
            unknown_sol = full(casadiFun.equilibrateMuscleTendonSingleMuscle(unknown_sol, known_sol)) ; % x find
            unknown_gast = full(casadiFun.equilibrateMuscleTendonSingleMuscle(unknown_gast, known_gast)) ; % x find
           
            % extract interest variales     
                % ta
            tendonForce_ta(compt) = unknown_ta(1) ;
            muscleForce_ta(compt) = unknown_ta(2) ;
            tendonLengthening_ta(compt) = unknown_ta(3) ;
            fiberLength_ta(compt) = unknown_ta(4)  ;
            pennationAngle_ta(compt) = unknown_ta(5)  ;
            tendonLength_ta(compt) = unknown_ta(3) + muscle_tendon_parameters_ta(4) ;
                %sol
            tendonForce_sol(compt) = unknown_sol(1) ;
            muscleForce_sol(compt) = unknown_sol(2) ;
            tendonLengthening_sol(compt) = unknown_sol(3) ;
            fiberLength_sol(compt) = unknown_sol(4)  ;
            pennationAngle_sol(compt) = unknown_sol(5)  ; 
            tendonLength_sol(compt) = unknown_sol(3) + muscle_tendon_parameters_sol(4) ;
                %gast
            tendonForce_gast(compt) = unknown_gast(1) ;
            muscleForce_gast(compt) = unknown_gast(2) ;
            tendonLengthening_gast(compt) = unknown_gast(3) ;
            fiberLength_gast(compt) = unknown_gast(4)  ;
            pennationAngle_gast(compt) = unknown_gast(5)  ;   
            tendonLength_gast(compt) = unknown_gast(3) + muscle_tendon_parameters_gast(4) ;

            % torque 
            torque(compt) = 1 ; % to do 
            
%                 if rand(1) > 0.9
%                         Force_length_representation(all_states_num,muscle_tendon_parameters_num,casadiFun,NameMuscles)
%                 end
        end
    end
end

            % output data 
%             Data = [torque'; a_num, ...
%                 length_UMT_ta(1,compt), length_UMT_sol(1,compt), length_UMT_gast(1,compt),...
%                 pennationAngle_ta, pennationAngle_sol, pennationAngle_gast,...
%                 tendonLength_ta, tendonLength_sol, tendonLength_gast] ; 
