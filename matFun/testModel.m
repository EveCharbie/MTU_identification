function testModel(skeleton_num,muscle_tendon_parameters_num,casadiFun)
Muscle_activation = .1 ; 
toleratedError = 0.1;
%% test of kinematics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % skeleton spatial configuration 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q1 = 0; % x
q2 = 0; % y
q3 = 0; % z
q4 = 0; % alpha hip
q5 = 0; % alpha knee
q6 = -25; % alpha ankle

% conversion deg to rad
q5 = (q5/180)*pi; q6 = (q6/180)*pi; 

q_num = [q1,q2,q3,q4,q5,q6];

    % plot skeleton 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Origin_num,Insertion_num,~,Markers_num] = casadiFun.ForwardKinematics([q_num,skeleton_num]);

plotmodel(Origin_num, Insertion_num, Markers_num)

    % mtu length
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mtu_length_num  = full(casadiFun.getUMTLength([q_num,skeleton_num]));

tibialis_length_num = mtu_length_num(1);
soleus_length_num = mtu_length_num(2);
gastrocnemius_length_num = mtu_length_num(3);

% disp(['tibials muscle-tendon unit length: ', num2str(tibialis_length_num),' m'])
% disp(['soleus muscle-tendon unit length: ', num2str(soleus_length_num),' m'])
% disp(['gastrocnemius muscle-tendon unit length: ', num2str(gastrocnemius_length_num),' m'])

    % mtu moment arm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mtu_moment_arm_num  = full(casadiFun.getMomentArm([q_num,skeleton_num]));

tibialis_moment_arm_num = mtu_moment_arm_num(1,end);
soleus_moment_arm_num = mtu_moment_arm_num(2,end);
gastrocnemius_moment_arm_num = mtu_moment_arm_num(3,end);
% 
% disp(['tibials muscle-tendon unit moment arm: ', num2str(tibialis_moment_arm_num),' m'])
% disp(['soleus muscle-tendon unit moment arm: ', num2str(soleus_moment_arm_num),' m'])
% disp(['gastrocnemius muscle-tendon unit moment arm: ', num2str(gastrocnemius_moment_arm_num),' m'])


    %% test rooted variables (signle muscle)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for tibialis
p_tibialis = muscle_tendon_parameters_num([1,4,7,10]);
[rooted_variables_tibialis,temp_residuals] = rootVariables(Muscle_activation,...
    tibialis_length_num,...
    p_tibialis,...
    casadiFun,...
    toleratedError,...
    1000);
pause

% for soleus
p_soleus = muscle_tendon_parameters_num([2,5,8,11]);
[rooted_variables_soleus,temp_residuals] = rootVariables(Muscle_activation,...
    soleus_length_num,...
    p_soleus,...
    casadiFun,...
    toleratedError,...
    1000);
pause

% for soleus
p_gastrocnemius = muscle_tendon_parameters_num([3,6,9,12]);
[rooted_variables_gastrocnemius,temp_residuals] = rootVariables(Muscle_activation,...
    gastrocnemius_length_num,...
    p_gastrocnemius,...
    casadiFun,...
    toleratedError,...
    1000);
pause

end 