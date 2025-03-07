function [Data,header] = HypotheticalDataGenerator(known_parameters_num,muscle_tendon_parameters_num,casadiFun)
%% 1. Set up the data generator 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output data name
header = {'Torque','q1','q2',...
    'activation_tibialis','activation_soleus','activation_gastrocnemius',...
    'fiber_tibialis','fiber_soleus','fiber_gastrocnemius',...
    'phi_tibialis','phi_soleus','phi_gastrocnemius'} ;

ntrialsFail = 0 ; % compteur of trials non-optimized (p > 10e-5)

ntrialsSucceds = 0 ; % compteur of trials optimized (p < 10e-5)

maxIteration = 100 ; % to find convergence 

% create a figure to observe the progressing 
fig =  uifigure("Name", "Data gerenrator", "Color", [1,1,1]) ;
fig.Position(3:4) = [550,150] ;
d = uiprogressdlg(fig,'Title','Please Wait',...
    'Message','Data generation');
d.Value = 0 ; 
pause(0.1)

%% 2. Selection of trials 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2.1 Musculo skeletical configuration during trial (input)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qankle = [-20 : 2.5 : 20 ]'...
    /180*pi ; % ankle angle 

qknee =  [0 : 5 : 85 , ...
    87.5 : -5 : 2.5]'...
    /180*pi ; % knee angle 

% 2.2 input : Neuronal activation (input)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_num = [0, 0, 0; ...
    0, 0.2, 0.2;...
    0, 0.4, 0.4;...
    0, 0.6, 0.6;...
    0.2, 0, 0;...
    0.4, 0, 0;...
    0.6, 0, 0]; % Muscle activation [Tibialis Anterior, Soleus, Gastrocnemius]

% 2.3 Muscle Tendon Parameters (ℓom, φo, Fom, ℓst)  (input)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
muscle_tendon_parameters_ta = muscle_tendon_parameters_num([1,4,7,10]) ; % Tibialis Anterior (ℓom, φo, Fom, ℓst)
muscle_tendon_parameters_sol = muscle_tendon_parameters_num([2,5,8,11]) ; % Soleus (ℓom, φo, Fom, ℓst)
muscle_tendon_parameters_gast = muscle_tendon_parameters_num([3,6,9,12]) ; % Gastrocnemius (ℓom, φo, Fom, ℓst)



%% 2. Data Generator 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trials = 0 ; % compt Iteration 
ntrials = size(a_num,1) * size(qknee,1) * size(qankle,1);   % total trials

for i = 1 : size(a_num,1)                                                  % for each muscle activation muscle
    for ii = 1 : size(qknee,1)                                             % for each knee angle
        for iii = 1 :  size(qankle,1)                                      % for each knee angle
            % 2.1 Progression
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            trials = trials + 1 ; % current trials 
            percent = round((trials/ntrials),2) ; 
            if d.Value < percent
                d.Value = percent ; 
                d.Message = ['Progressing : ', num2str(percent*100),' %'];
                pause(0.1)
            end

            % 2.2 initial guess
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            unknown_ta = [muscle_tendon_parameters_ta(4) ,muscle_tendon_parameters_ta(1:2)] ;
            unknown_sol = [muscle_tendon_parameters_sol(4) ,muscle_tendon_parameters_sol(1:2)] ;
            unknown_gast = [muscle_tendon_parameters_gast(4) ,muscle_tendon_parameters_gast(1:2)] ;

            % 2.3 Neuromusculoskeletal configuration
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            q_num = [0, 0 , 0, 0, qknee(ii), qankle(iii) ] ; 
            musculoskeletal_states_num = [q_num , known_parameters_num] ;

            % 2.4 UMT length
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            UMT_Length = full(casadiFun.getUMTLength(musculoskeletal_states_num)) ; 

            % 2.5 Solver 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 2.5.1 Tibialis Anterior optimisation 
            known_ta = [a_num(i,1), UMT_Length(1) , muscle_tendon_parameters_ta] ; % F = f(a, Lmtu,p[ℓom, φo, Fom, ℓst])
            % unknown_sol(1) = UMT_Length(1)- (muscle_tendon_parameters_ta(1))*cos(muscle_tendon_parameters_ta(2))
            unknown_ta = full(casadiFun.equilibrateMuscleTendonSingleMuscle2(unknown_ta, known_ta)) ; % equilibrium            
            err_ta = casadiFun.equilibriumErrorSingleMuscle2(unknown_ta',known_ta); % residuals
            
            if any(full(err_ta)>1e-5) % if too much residuals
                [unknown_ta,err_ta] = NewStart(casadiFun,known_ta,unknown_ta,err_ta,maxIteration,1e-5) ; 
            end

            % 2.5.2 Soleus  optimisation
            known_sol = [a_num(i,2), UMT_Length(2) , muscle_tendon_parameters_sol] ;
            unknown_sol = full(casadiFun.equilibrateMuscleTendonSingleMuscle2(unknown_sol, known_sol)) ; % x equilibrium
            err_sol = casadiFun.equilibriumErrorSingleMuscle2(unknown_sol',known_sol);

            if any(full(err_sol)>1e-5) % if too much residuals
                [unknown_sol,err_sol] = NewStart(casadiFun,known_sol,unknown_sol,err_sol,maxIteration,1e-5) ;
            end

            % 2.5.3 Gastrocnemius anterior optimisation
            known_gast = [a_num(i,3), UMT_Length(3) , muscle_tendon_parameters_gast] ;
            unknown_gast = full(casadiFun.equilibrateMuscleTendonSingleMuscle2(unknown_gast, known_gast)) ; % x equilibrium
            err_gast = casadiFun.equilibriumErrorSingleMuscle2(unknown_gast',known_gast);

            if any(full(err_gast)>1e-5) % if too much residuals
                [unknown_gast,err_gast] = NewStart(casadiFun,known_gast,unknown_gast,err_gast,maxIteration,1e-5) ;
            end

            % 2.6 Check if there is too mutch error
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if any([full(err_ta),full(err_sol),full(err_gast)] >1e-5) % if too much residuals
                ntrialsFail = ntrialsFail + 1 ;
                continue
            else
                ntrialsSucceds = ntrialsSucceds + 1 ; 
            end

            % 2.7 Data extraction (to improve)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Musculo skeletical configuration
            q1(ntrialsSucceds) = qknee(ii) ;
            q2(ntrialsSucceds) = qankle(iii) ;

            % Neuronal activation
            a_ta(ntrialsSucceds) = a_num(i,1);
            a_sol(ntrialsSucceds) = a_num(i,2);
            a_gast(ntrialsSucceds) = a_num(i,3);

            % Muscle Tendon Unit length 
            length_UMT_ta(1,ntrialsSucceds) = UMT_Length(1);
            length_UMT_sol(1,ntrialsSucceds) = UMT_Length(2) ;
            length_UMT_gast(1,ntrialsSucceds) = UMT_Length(3) ;

            % rooted variables 
                % Tibialis Anterior
            tendonLength_ta(ntrialsSucceds) = unknown_ta(1) ;
            fiberLength_ta(ntrialsSucceds) = unknown_ta(2)  ;
            pennationAngle_ta(ntrialsSucceds) = unknown_ta(3)  ;

                % Soleus
            tendonLength_sol(ntrialsSucceds) = unknown_sol(1) ;
            fiberLength_sol(ntrialsSucceds) = unknown_sol(2)  ;
            pennationAngle_sol(ntrialsSucceds) = unknown_sol(3)  ; 

                % Gastrocnemius
            tendonLength_gast(ntrialsSucceds) = unknown_gast(1) ;
            fiberLength_gast(ntrialsSucceds) = unknown_gast(2)  ;
            pennationAngle_gast(ntrialsSucceds) = unknown_gast(3)  ; 

            %% def
            % rooted
            %%%%%%%%%%%
            % fiberLength [3]
            % tendonLength [3]
            % rootedvariables = [fiberLength, tendonLength] [6]

            % concatenated
            %%%%%%%%%%%
            % "musculoskeletal_states" = [q, known_parameter] [39]
            % "neuromusculoskeletal_state" = [a, q, known_parameter] [42]
            % "all_states" = [neuromusculoskeletal_state, rooted_variales] [48]

            % about forces
            neuromusculoskeletal_states = [a_ta(ntrialsSucceds),a_sol(ntrialsSucceds),a_gast(ntrialsSucceds),musculoskeletal_states_num];

            rooted_variables = [fiberLength_ta(ntrialsSucceds),fiberLength_sol(ntrialsSucceds),fiberLength_gast(ntrialsSucceds),...
                tendonLength_ta(ntrialsSucceds),tendonLength_sol(ntrialsSucceds),tendonLength_gast(ntrialsSucceds)];

            all_states_num = [neuromusculoskeletal_states,rooted_variables];

            torque(ntrialsSucceds) = full(casadiFun.getJointMoment(all_states_num,muscle_tendon_parameters_num));
        end
    end
end
close(fig)

%% Output of the function 
% output data
pennationAngle_ta = mod(pennationAngle_ta,pi);
pennationAngle_sol = mod(pennationAngle_sol,pi);
pennationAngle_gast = mod(pennationAngle_gast,pi);

Data = [torque', q1', q2',... 1..3
    a_ta', a_sol', a_gast' ... 4..6
    length_UMT_ta', length_UMT_sol', length_UMT_gast',... 7..9
    fiberLength_ta', fiberLength_sol', fiberLength_gast',... 10..12
    pennationAngle_ta', pennationAngle_sol', pennationAngle_gast'] ; %25..27

% save data 
save("Data.mat", "Data")

% display in command Windox the score 
disp('fail trials in percent :')
%disp(round((ntrialsFail/ntrials)*100,2))
disp((ntrialsFail/ntrials)*100)

end


function [unknown,err] = NewStart(casadiFun,known,unknown,err,maxIteration,objective)
% when the solver fails to converge. The solution that did not converge properly
temp = unknown;
comptErr = 0 ;
while comptErr < maxIteration
    random_values = 0.2 * randn(1, size(unknown,1)) + 1;               % Generate random values from a normal distribution
    % random_values(random_values < 0.5) = 0.5; random_values(random_values > 1.5) = 1.5; % between 50 % and 150%
    random_values(random_values < 0.5) = 0.8; random_values(random_values > 1.5) = 1.2; % between 50 % and 150%

    unknown = unknown .* random_values' ;
    unknown = full(casadiFun.equilibrateMuscleTendonSingleMuscle2(unknown, known)) ; % equilibrium
    err = casadiFun.equilibriumErrorSingleMuscle2(unknown',known);  % residuals

    if all(full(err)<objective)
        break
    end
    comptErr = comptErr + 1 ;
end

% if it dit not converge 
if all(full(err)>objective)
    unknown = temp;
end
end 