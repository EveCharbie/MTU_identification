%% test casadi functions
%run('.\De_Groote'); 

example_state_param; 

[muscle_origin_in0, muscle_insertion_in0, joint_centers] = ForwardKinematics(q_num, known_parameters_num);
umt_length = UMTLength(q_num, known_parameters_num);
moment_arms = MomentArm(q_num, known_parameters_num);


% plotmodel(muscle_origin_in0, muscle_insertion_in0, joint_centers)


%% test du model réologique

lm_num = [0 : 0.01 : 1.7 ; 0 : 0.01 : 1.7 ; 0 : 0.01 : 1.7] ; 
lt_num =  [linspace(0.95,1.05,length(lm_num(1,:))) ; linspace(0.95,1.05,length(lm_num(1,:))) ; linspace(0.95,1.05,length(lm_num(1,:)))] ; 
a_num = ones(3,1) ;

% from OISM
% % soleus
% Sol_f0m = 3600 ;
% Sol_l0m = 0.031 ;
% Sol_phi0 = 0.20943951 ; %rad
% Sol_lst = 0.31 ;
% 
% % Gast
% Gast_f0m = 2500 ;
% Gast_l0m = 0.09 ;
% Gast_phi0 = 0.29670597 ; %rad
% Gast_lst = 0.36 ;
% 
% % Ta
% Ta_f0m = 3000 ;
% Ta_l0m = 0.098 ;
% Ta_phi0 = 0.29670597 ; %rad
% Ta_lst = 0.08726646 ;
% 
% % 2. Déterminer des bornes physiologiques pour les paramètres (en faire des vecteurs : min_param et max_param)
% %min_param =
% %Ta_phi0 = 0
% % Ta_l0m = 0
% %max_param =
% % Ta_phi0 = pi/2
% % Ta_l0m = 0.25
% % Ta_lst =
%%
for i = 1 :  length (lm_num)
    %% Normalized 
    % passives muscle force (Normalized)
    PMFN = Normalized_Muscle_Passive_Force_Length(q_num,known_parameters_num,lm_num(:,i)) ; 
    NormalizedMusclePassiveForceLength_num(:,i) = full(PMFN) ; 
    % active muscle force (Normalized)
    AMFN = Normalized_Muscle_Active_Force_Length(q_num,known_parameters_num,lm_num(:,i)) ; 
    NormalizedMuscleActiveForceLength_num(:,i) = full(AMFN) ; 
    % total muscle force (Normalized)
    MFN = Normalized_Muscle_Force(q_num,known_parameters_num,lm_num(:,i),a_num(:)) ; 
    NormalizedMuscleForce_num(:,i) = full(MFN) ; 
    % total tendon force (Normalized)
    TFN = Normalized_Tendon_Force(q_num,known_parameters_num,lt_num(:,i)) ; 
    NormalizedTendonForce_num(:,i) = full(TFN) ; 

    %% Non-Normalized
    % passives muscle force 
    PMF = Muscle_Passive_Force_Length(q_num,known_parameters_num,muscle_tendon_parameters_num,lm_num(:,i)) ;
    MusclePassiveForceLength_num(:,i) = full(PMF) ;
    % active muscle force
    AMF = Muscle_Active_Force_Length(q_num,known_parameters_num,muscle_tendon_parameters_num,lm_num(:,i),a_num(:)) ;
    MuscleActiveForceLength_num(:,i) = full(AMF) ;
    % total muscle force 
    MF = Muscle_Force(q_num,known_parameters_num,muscle_tendon_parameters_num,lm_num(:,i),a_num(:)) ;
    MuscleForce_num(:,i) = full(MF) ;
    % total tendon force 
    TF = Tendon_Force(q_num,known_parameters_num,muscle_tendon_parameters_num,lt_num(:,i)) ;
    TendonForce_num(:,i) = full(TF) ;

end 

%% illustration
%% relation force longueur Normalized

figure("Name","relation force longueur ","Color",[1 1 1])
subplot(1,2,1)
plot(lm_num(1,:) , NormalizedMusclePassiveForceLength_num(1,:), 'k')
hold on 
plot(lm_num(1,:) , NormalizedMuscleActiveForceLength_num(1,:), 'r')
plot(lm_num(1,:) , NormalizedMuscleForce_num(1,:), 'k.')
xlabel('Normalized length','FontWeight','bold')
ylabel('Normalized force','FontWeight','bold')
xlim([0.2 1.7])
ylim([0 1.4])
title('Muscle','FontAngle','italic')
legend(["Passive", "Active","Total"])
hold off

subplot(1,2,2)
plot(lt_num(1,:) , NormalizedTendonForce_num(1,:), 'k')
hold on 
xlabel('Normalized length','FontWeight','bold')
ylabel('Normalized force','FontWeight','bold')
xlim([0.95 1.05])
ylim([0 1])

title('Tendon','FontAngle','italic')
hold off


%% relation force longueur 
figure("Name","relation force longueur ","Color",[1 1 1])
subplot(1,2,1)
plot(lm_num(1,:) , MusclePassiveForceLength_num(1,:), 'k')
hold on 
plot(lm_num(1,:) , MuscleActiveForceLength_num(1,:), 'r')
plot(lm_num(1,:) , MuscleForce_num(1,:), 'k.')
xlabel('Normalized length','FontWeight','bold')
ylabel('Force (N)','FontWeight','bold')
xlim([0.2 1.7])
%ylim([0 1.4])
title('Muscle','FontAngle','italic')
legend(["Passive", "Active","Total"])
hold off

subplot(1,2,2)
plot(lt_num(1,:) , TendonForce_num(1,:), 'k')
hold on 
xlabel('Normalized length','FontWeight','bold')
ylabel('Force (N)','FontWeight','bold')
xlim([0.95 1.05])
%ylim([0 1])

title('Tendon','FontAngle','italic')
hold off




%% test of the  moment arm function
for i = 1 :  length (lm_num)
    MA  = Momentarticualire(q_num,known_parameters_num,muscle_tendon_parameters_num,lt_num(:,i)); 
    Moment(:,i) = full(MA) ; 
end

figure
plot(lt_num(1,:),Moment(5,:))


