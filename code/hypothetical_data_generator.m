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
qankle = linspace(-30,30, 10)/180*pi ;
qknee = [0, 10, 20, 30, 40, 50, 60, 70, 80]/180*pi ;

% input :  muscle activation  (random)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_num = [0, 0, 0; ...
    0, 1, 1;...
    1, 0, 0;];

ntrials = 0 ;
ntrialsfail = 0 ; 
Data = [];

%%

unknown_num = [ones(1,9).* .0001 ,muscle_tendon_parameters_num(1:6) ] ;

    compt = 0 ; 
for i = 1 : size(a_num,1) % activation muscle
    for ii = 1 : length(qknee)
        for iii = 1 :  length(qankle)
            q_num = [0, 0 , 0, 0, qknee(ii), qankle(iii) ] ; %
            musculoskeletal_states_num = [q_num , known_parameters_num] ;
            neuromusculoskeletal_states_num = [a_num(i,:), musculoskeletal_states_num] ; 

            p_num = horzcat( neuromusculoskeletal_states_num, muscle_tendon_parameters_num) ;
            
            compt = compt+1 ; 
            temp = casadiFun.getUMTLength(musculoskeletal_states_num) ; 
            lta(1,compt) = full( temp(1)) ;
            lsol(1,compt) = full( temp(2)) ;
            lgast(1,compt) = full( temp(3)) ;
            

            known_num = [a_num(i,:), full(temp)' , muscle_tendon_parameters_num] ; 
            x_num = full(casadiFun.equilibrateMuscleTendon1(unknown_num, known_num)) ; % x find
            
            try
                x_num = full(casadiFun.equilibrateMuscleTendon1(unknown_num, known_num)) ; % x find

            catch
                 x_num = nan(15,1) ; 
            end

            if any(isnan(x_num)) %any(isnan(x_num))
                % do nothing 
                ntrialsfail = ntrialsfail +1 ; 
                idx(ntrialsfail,:) = [i, ii, iii] ; 

            else
                tendonForce_num = x_num(1:3)' ;
                muscleForce_num = x_num(4:6)' ;
                tendonLengthening_num = x_num(7:9)' ;
                fiberLength_num = x_num(10:12)' ;
                pennationAngle = x_num(13:15)' ;

                rootedvariables =  [fiberLength_num,tendonLengthening_num] ;
                all_states_num = [neuromusculoskeletal_states_num, rootedvariables] ;
                Torque = full(casadiFun.getJointMoment(all_states_num,muscle_tendon_parameters_num)) ;                                             % torque
                
                % output
                ntrials = ntrials +1 ;
                Data(ntrials,:) = [Torque, qknee(ii) , qankle(iii),  a_num(i,:),...
                    fiberLength_num, pennationAngle, tendonLengthening_num ] ; % variables mesured

                if i == 2
% %                    1+1
% %                     [Origin,Insertion,ViaPoint,Markers]  = casadiFun.ForwardKinematics(musculoskeletal_states_num) ; 
% %                     plotmodel(Origin, Insertion, Markers)
% %                     Force_length_representation(all_states_num,muscle_tendon_parameters_num,casadiFun,NameMuscles)
                end 
%                 if rand(1) > 0.9
%                         Force_length_representation(all_states_num,muscle_tendon_parameters_num,casadiFun,NameMuscles)
%                 end
            end
        end
    end
end
%%
% color = [linspace(0.1,0.9,length(qknee)); ...
%     linspace(0.1,0.9,length(qknee));...
%     linspace(0.1,0.9,length(qknee))] ; 
% 
% figure
% hold on 
% lgd = {} ; 
% 
% for iiii = 1 : length(qknee)
% 
% plot(((qankle/pi)*180),lgast(iiii,:),'DisplayName','lgast', 'Color',color(:,iiii)')
% lgd{iiii} = ['q knee : ',num2str((qknee(iiii)/pi)*180)] ;
% end
% legend(lgd)
% xlabel('q ankle (°)','FontWeight','bold')
% ylabel('UMT length (m)','FontWeight','bold')
% 
% fprintf('fail root  :\n')
% fprintf(num2str(ntrialsfail))
% fprintf('\n')
% fprintf('fail percent  :\n')
% fprintf(num2str((ntrialsfail/ntrials) *100 ))
% fprintf('\n')