%% test casadi functions
%run('.\De_Groote'); 
%% test du model géometrique (OK)
De_Groote ; 
% foot; leg; thigh
segment_length_num = [.2; .4; .45];

% Origin TA, SO, GA; insertion TA; SO; GA
muscle_origin = [-.3; .01 ; 0 ; 
    -.3 ; -.01 ; 0 ; 
    -.02; -.01 ; 0 ];

muscle_insersion = [-.17; 0.01; 0 ; 
    -.23; -.02; 0 ; 
    -.23; -.02; 0];

known_parameters_num = [segment_length_num;muscle_origin;muscle_insersion];

q_num = [0; 0; 0; -20/180*pi; 10/180*pi; 50/180*pi];


[muscle_origin_in0, muscle_insertion_in0, joint_centers] = ForwardKinematics(q_num, known_parameters_num);
umt_length = UMTLength(q_num, known_parameters_num);
moment_arms = MomentArm(q_num, known_parameters_num);


% plotmodel(muscle_origin_in0, muscle_insertion_in0, joint_centers)


%% test du model réologique

lm_num = [0 : 0.01 : 1.7 ; 0 : 0.01 : 1.7 ; 0 : 0.01 : 1.7] ; 
lt_num =  [linspace(0.95,1.05,length(lm_num(1,:))) ; linspace(0.95,1.05,length(lm_num(1,:))) ; linspace(0.95,1.05,length(lm_num(1,:)))] ; 
a_num = ones(3,1) ;
% TA; SO; GA;
l0m_num = [0.098 ;  0.031 ; 0.09] ; 
phi0_num = [0.29670597 ; 0.20943951 ; 0.29670597 ] ;
f0m_num = [3000 ; 3600 ; 2500] ; 
lst_num = [0.08726646 ; 0.31 ; 0.36 ] ; 
muscle_tendon_parameters_num =  [l0m_num ; phi0_num ; f0m_num ; lst_num ] ; 

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


%% faire longueur non nomralized 
%% fonction de contrainte 
length_tendon = lst_num + 0.01;
length_muscle = SX.sym('Muscle_length',3) ;

length_tendon_normalized_num = Length_Tendon_Nomralized(length_tendon,muscle_tendon_parameters_num) ; 
length_muscle_normalized_num = Length_Muscle_Nomralized(length_muscle,muscle_tendon_parameters_num) ; 

Ttest = Tendon_Force(q_num,known_parameters_num,muscle_tendon_parameters_num,length_tendon_normalized_num) ;


g0 = TendonForce - (cos(pennation_angle) .* MuscleForce) ; 
g1 = umt_length' - (cos(pennation_angle) .* length_muscle + length_tendon) ; 

g = Function('g',[length_tendon,length_muscle],[g0,g1]) ; 
G = rootfinder('G','newton',g) ; 

%%
% z = SX.sym('x',nz);
% x = SX.sym('x',nx);
% 
% g = Function('g',{z,x},{g0,g1});
% G = rootfinder('G','newton',g);
