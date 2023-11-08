%% hypothetical data
nTrials = 6 ;

% passives triales systeme configuration (muscle activation and q's)

rang_ankle = linspace(-40,50, 1000) ; % range ankle
qankle = ones(nTrials, length(rang_ankle)) .* rang_ankle;
a_ = zeros(6,3);
% a_ = ones(6,3) .* .2 ; 
% a_(:,1) = zeros(6,1) ; 
% a_(:,2:3) = zeros(6,2) ;
qknee  = [0, 20, 40, 60, 80, 100];


for i = 1 : nTrials
    for ii = 1 :  length(rang_ankle)
        %% input of root finder function
        q_num = [0, 0 , 0, 0, (qknee(i)/180)*pi , (qankle(i,ii)/180)*pi ] ;

        [muscle_origin_in0, muscle_insertion_in0, joint_centers] = ForwardKinematics(q_num, known_parameters_num);
        umt_length = getUMTLength(q_num, known_parameters_num) ;
        moment_arms = getMomentArm(q_num, known_parameters_num) ;
        musculoskeletal_states_num = [q_num , known_parameters_num] ;

        p_num = horzcat(a_(i,:),musculoskeletal_states_num , muscle_tendon_parameters_num) ;

        %% root finder
        x0 = [0.001, 0.001, 0.001 ,l0m_num] ;  % tendonLengthening, fibre length   (TA SOL GAST)

        x_num = equilibrateMuscleTendon(x0, p_num) ;
        try
            tendonLengthening_num = x_num(1:nMuscles)';
            fiberLength_num = x_num(nMuscles+1:end)';

            all_states_num = [musculoskeletal_states_num,  fiberLength_num, tendonLengthening_num , a_num] ;
            Mtemp = full(getJointMoment(all_states_num,muscle_tendon_parameters_num)  );
            FMtemp = full(getMusclePassiveForce(all_states_num,muscle_tendon_parameters_num)  );
            Moment(i,ii) = Mtemp(end) ; 
            FMtemp(i,ii) = FMtemp(end) ; 
        catch
            Moment(i,ii) = nan(1) ;
            FMtemp(i,ii) = nan(1) ;
        end
    end

end
figure 
plot(rang_ankle,Moment')
text(-15, min(min(Moment)) - 0.1, 'dorsiflexion', 'HorizontalAlignment', 'center');
text(25,  min(min(Moment)) , 'plantar flexion', 'HorizontalAlignment', 'center');
% figure 
% plot(rang_ankle,FMtemp')