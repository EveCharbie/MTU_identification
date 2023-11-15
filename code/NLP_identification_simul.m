
% Call data_generator




% Add noise


% Start with an empty NLP
w={}; %variables
w0 = []; % initial guess
lbw = []; % lower bound of variable
ubw = []; % upper bound of variable
J = 0; % cost | objective function
g={}; % constraints
lbg = []; % lower bound of constraints
ubg = []; % upper bound of constraints


% create muscle parameters variables
w = { w{:}, unknown_parameters};

% we should find muscle_tendon_parameters_num
w0 =  [w0; muscle_tendon_parameters_num + rand()*xxx]; %todo

lbw = [lbw; muscle_tendon_parameters_num * 0.5];%todo
ubw = [ubw; muscle_tendon_parameters_num * 2];%todo



% header = {'Torque','q1','q2',... 1..3
%     'activation_tibialis','activation_soleus','activation_gastrocnemius',...
%     'momentarm_tibialis','momentarm_soleus','momentarm_gastrocnemius',...
%     'UMT_tibialis','UMT_soleus','UMT_gastrocnemius',... todo: REMOVE
%     'fiber_tibialis','fiber_soleus','fiber_gastrocnemius',...
%     'phi_tibialis','phi_soleus','phi_gastrocnemius',...
%      'tendon_tibialis','tendon_soleus','tendon_gastrocnemius'} ;


% create variable for each trial

% Weightings in cost function
W_torque = 1;
W_length = xx; %todo 1Nm correspond à 1 mm?
W_angle = xx;  % 1 Nm correspond à 3 deg?

for trial = 1:size(results, 1)
    fiberLength_k = MX.sym(['Fiber_length_' num2str(trial)], nMuscles);
    tendonLengthening_k = MX.sym(['Tendon_Lengthening_' num2str(trial)], nMuscles);

    w = { w{:}, fiberLength_k, tendonLengthening_k};
    w0 =  [w0; xxxx]; %todo add initial guess
    lbw = [lbw; ones(3,1)*0.01; zeros(3,1)];%todo
    ubw = [ubw; ones(3,1)*0.4; ones(3,1)*.2];%todo


    data = trial{trial};

    % all_states = vertcat(musculoskeletal_states, muscle_tendon_states);
    % musculoskeletal_states = vertcat(q, known_parameters)
    % muscle_tendon_states = vertcat(fiberLength, tendonLengthening, a);

    q = data(2:3); 

    all_states_num = [q; known_parameters; musculoskeletal_states; muscle_tendon_states]; %todo
    
    


    
    
    umt_lengths = getUMTLength(q, known_parameters);
    torque = getJointMoment(xxx); %todo add input
    tendonL = 
    
    % constraints
    contraints = g(x,p); %todo rename g in de_groot
    g={g, constraints}; 
    lbg = [lbg; zeros(6,1)]; 
    ubg = [ubg; zeros(6,1)]; 

    % objective


    J = J + W_torque * (data(1) - torque)^2; %add error on joint torque
    
    % add error on tendon length and pennation angle
    J = J + W_length * (data(xxx) - tendonLength).^2;
    J = J + W_angle * (data(xxx) - pennationAngle).^2;

end






% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob);

% Solve the NLP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
            'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);

param_opt = w_opt(1:nparam);

%compare with muscle_tendon_parameters_num
