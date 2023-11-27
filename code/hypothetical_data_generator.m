%% hypothetical data
% import model

clear 
close all
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

% input : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_num = [0, 0, 0; ...
    0, 0.2, 0.2; ...
    0, 0.4, 0.4;...
    0, 0.6, 0.6;...
    0, 0.8, 0.8;...
    0, 1, 1;...
    0, 0.9, 0.9;...
    0, 0.7, 0.7;...
    0, 0.5, 0.5;...
    0, 0.3, 0.3;...
    0, 0.1, 0.1;...
    0.2, 0, 0; ...
    0.4, 0, 0;...
    0.6, 0, 0;...
    0.8, 0, 0;...
    1, 0, 0];

ntrials = 0 ;
ntrialsfail = 0 ; 
Data = [];

    % muscle tendon parameter 
muscle_tendon_parameters_ta = muscle_tendon_parameters_num([1,4,7,10]) ; 
muscle_tendon_parameters_sol = muscle_tendon_parameters_num([2,5,8,11]) ; 
muscle_tendon_parameters_gast = muscle_tendon_parameters_num([3,6,9,12]) ; 
    % inital guess
unknown_ta = [ones(1,3).* .0001 ,muscle_tendon_parameters_ta(1:2)] ;
unknown_sol = [ones(1,3).* .0001 ,muscle_tendon_parameters_sol(1:2)] ;
unknown_gast = [ones(1,3).* .0001 ,muscle_tendon_parameters_gast(1:2)] ;


    compt = 0 ; 
    compt2 = 0;
for i = 1 : size(a_num,1) % activation muscle
    for ii = 1 : length(qknee)
        for iii = 1 :  length(qankle)

            compt = compt + 1;
            % Neuromusculoskeletal configuration
            q_num = [0, 0 , 0, 0, qknee(ii), qankle(iii) ] ; 
            musculoskeletal_states_num = [q_num , known_parameters_num] ;
            neuromusculoskeletal_states_num = [a_num(i,:), musculoskeletal_states_num] ; 
            p_num = horzcat( neuromusculoskeletal_states_num, muscle_tendon_parameters_num) ;

           %input of the equilibrium function 

                % UMT length 
            UMT_Length = full(casadiFun.getUMTLength(musculoskeletal_states_num)) ; 

            
                % known parameter 
            known_ta = [a_num(i,1), UMT_Length(1) , muscle_tendon_parameters_ta] ;
            known_sol = [a_num(i,2), UMT_Length(2) , muscle_tendon_parameters_sol] ;
            known_gast = [a_num(i,3), UMT_Length(3) , muscle_tendon_parameters_gast] ;

            % Equilibrium
            unknown_ta = full(casadiFun.equilibrateMuscleTendonSingleMuscle(unknown_ta, known_ta)) ; % x find            
            unknown_sol = full(casadiFun.equilibrateMuscleTendonSingleMuscle(unknown_sol, known_sol)) ; % x find
            unknown_gast = full(casadiFun.equilibrateMuscleTendonSingleMuscle(unknown_gast, known_gast)) ; % x find
    

            %test
            x = [unknown_ta, unknown_sol, unknown_gast];
            x=x'; x=x(:);

            y = [known_ta; known_sol; known_gast];
            y=y(:);
            err = casadiFun.equilibriumError1(x,y);
            
            % random initial guess
            if any(full(err)>1e-5)
                compterr = 0 ;
                while compterr < 30
                    xtemp = x * rand(1) ; 


                    
                    err = casadiFun.equilibriumError1(xtemp,y);

                    if any(full(err)<1e-5)
                        x = xtemp ; 
                        break
                    end
                    compterr = compterr + 1 ;
                end
                if any(full(err)>1e-5)
                    continue
                end
                fprintf('ERROR \n')
            end

            compt2 = compt2+1 ; 



            % extract interest variales     
                % ta
            q1(compt2) = qknee(ii) ; 
            q2(compt2) = qankle(iii) ;


            a_ta(compt2) = a_num(i,1);
            a_sol(compt2) = a_num(i,2);
            a_gast(compt2) = a_num(i,3);


            length_UMT_ta(1,compt2) = UMT_Length(1);
            length_UMT_sol(1,compt2) = UMT_Length(2) ;
            length_UMT_gast(1,compt2) = UMT_Length(3) ;            

            tendonForce_ta(compt2) = unknown_ta(1) ;
            muscleForce_ta(compt2) = unknown_ta(2) ;
            tendonLengthening_ta(compt2) = unknown_ta(3) ;
            fiberLength_ta(compt2) = unknown_ta(4)  ;
            pennationAngle_ta(compt2) = unknown_ta(5)  ;
            tendonLength_ta(compt2) = unknown_ta(3) + muscle_tendon_parameters_ta(4) ;
                %sol
            tendonForce_sol(compt2) = unknown_sol(1) ;
            muscleForce_sol(compt2) = unknown_sol(2) ;
            tendonLengthening_sol(compt2) = unknown_sol(3) ;
            fiberLength_sol(compt2) = unknown_sol(4)  ;
            pennationAngle_sol(compt2) = unknown_sol(5)  ; 
            tendonLength_sol(compt2) = unknown_sol(3) + muscle_tendon_parameters_sol(4) ;
                %gast
            tendonForce_gast(compt2) = unknown_gast(1) ;
            muscleForce_gast(compt2) = unknown_gast(2) ;
            tendonLengthening_gast(compt2) = unknown_gast(3) ;
            fiberLength_gast(compt2) = unknown_gast(4)  ;
            pennationAngle_gast(compt2) = unknown_gast(5)  ;   
            tendonLength_gast(compt2) = unknown_gast(3) + muscle_tendon_parameters_gast(4) ;

            % torque 
            torque(compt2) = full(casadiFun.getJointMoment2(musculoskeletal_states_num, ...
                [tendonForce_ta(compt2),tendonForce_sol(compt2),tendonForce_gast(compt2)] ) ); 
            
%                 if rand(1) > 0.9
%                         Force_length_representation(all_states_num,muscle_tendon_parameters_num,casadiFun,NameMuscles)
%                 end
        end
    end
end

% output data
Data = [torque', q1', q2',... 1..3
    a_ta', a_sol', a_gast' ... 4..6
    length_UMT_ta', length_UMT_sol', length_UMT_gast',... 7..9
    fiberLength_ta', fiberLength_sol', fiberLength_gast',... 10..12
    pennationAngle_ta', pennationAngle_sol', pennationAngle_gast',... 13..15
    tendonLength_ta', tendonLength_sol', tendonLength_gast'... 16..18
    tendonForce_ta', tendonForce_sol', tendonForce_gast',... 19..21
    muscleForce_ta', muscleForce_sol', muscleForce_gast', ...22..24
    tendonLengthening_ta', tendonLengthening_sol', tendonLengthening_gast'] ; %25..27

% Best start : a   qknee   qankle  FT  FM  tendonLengthening   fiberLength pennationAngle
BestStart_ta = [a_ta', length_UMT_ta', tendonForce_ta', muscleForce_ta',...
    tendonLengthening_ta', fiberLength_ta', pennationAngle_ta'] ; 

BestStart_sol = [a_sol',length_UMT_sol', tendonForce_sol', muscleForce_sol',...
    tendonLengthening_sol', fiberLength_sol', pennationAngle_sol'] ; 

BestStart_gast = [a_gast', length_UMT_gast', tendonForce_gast', muscleForce_gast',...
    tendonLengthening_gast', fiberLength_gast', pennationAngle_gast'] ; 


save("Data.mat", "Data")
