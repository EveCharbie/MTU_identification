% import casadi
[ret, name] = system('hostname');
if strcmp(name, '942-27984')
    addpath('C:\Users\Stage\Desktop\Neuromusculoskeletal Modeling\Casadi')
elseif strcmp(name(1:27), 'MacBook-Air-de-mickaelbegon')
    addpath('/Users/mickaelbegon/Downloads/casadi-3.6.3-osx64-matlab2018b/')
end

import casadi.*


% Call data_generator
hypothetical_data_generator ;

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

UP = vertcat(unknown_parameters(:,1),unknown_parameters(:,2),...
    unknown_parameters(:,3),unknown_parameters(:,4));

% create muscle parameters variables
w = { w{:}, UP};
% we should find muscle_tendon_parameters_num
w0 =  [w0; muscle_tendon_parameters_num' ]; %TODO: add noise
lbw = [lbw; muscle_tendon_parameters_num' * 0.5]; % lower bound of variable
ubw = [ubw; muscle_tendon_parameters_num' * 2]; % upper bound of variable


% header = header = {'Torque','q1','q2',...
%     'activation_tibialis','activation_soleus','activation_gastrocnemius',...
%     'fiber_tibialis','fiber_soleus','fiber_gastrocnemius',...
%     'phi_tibialis','phi_soleus','phi_gastrocnemius'} ;


% create variable for each trial

% Weightings in cost function
W_torque = 1; % 1 Nm
W_length = 1;% / 0.005; %todo 1Nm correspond à 5 mm
W_angle = 1;% / ((3/180) * pi);  % 1 Nm correspond à 3 deg

eTorque = [];
eFiber = [];
ePennation = [];

for trial = 1:ntrials % for 1 to nb trials
    data = Data(trial,:); % known variables 

    fiberLength_k = SX.sym(['Fiber_length_' num2str(trial)], nMuscles);
    tendonLengthening_k = SX.sym(['Tendon_Lengthening_' num2str(trial)], nMuscles);
    w = { w{:}, tendonLengthening_k, fiberLength_k};
    w0 =  [w0; data(end-2:end)'; data(7:9)']; 
%     w0 =  [w0; ones(3,1)*0.001; ones(3,1)*0.1]; 
    lbw = [lbw; zeros(3,1); ones(3,1)*0.01]; % lower bound of variable
    ubw = [ubw; ones(3,1)*0.5; ones(3,1)*.07];    


    a_trial = data(4:6) ; % muscle activation during the trial
    q_trial = [0, 0, 0, 0, data(2:3)] ; % skeleton configuration during the trial
    musculoskeletal_states_trial = [q_trial, known_parameters_num ] ; % muscleskeleton configuration during the trial
    neuromusculoskeletal_states_trial = [a_trial, musculoskeletal_states_trial] ; % Neuromusculoskeletal states
    p_trial = vertcat(neuromusculoskeletal_states_trial', UP) ;


    %OK: UP, equilibriumError, 
    
    % constraints
    constraints = casadiFun.equilibriumError( ...
        vertcat(tendonLengthening_k, fiberLength_k), ...
        p_trial);

    g = { g{:}, constraints};
    lbg = [lbg; zeros(6,1)]; % lower bound of constraints
    ubg = [ubg; zeros(6,1)]; % upper bound of constraints
    
    all_states_trial = vertcat( ...
        neuromusculoskeletal_states_trial', ...
        tendonLengthening_k, ...
        fiberLength_k ) ;

    temp = casadiFun.getJointMoment(all_states_trial,  UP) ;
    Torque_simulated = temp(end) ;
    phi_simulated = casadiFun.getPennationAngle(all_states_trial, UP)' ;

    % objective
    eTorque = [eTorque; data(1) - Torque_simulated];
    eFiber = [eFiber; data(7:9)' - fiberLength_k];
    ePennation = [ePennation; data(10:12) - phi_simulated];

    J = J + W_torque * eTorque(trial)^2; %add error on joint torque
    J = J + W_length * sum(eFiber(trial,:).^2);% add error on tendon length and pennation angle
    J = J + W_angle * sum(ePennation(trial,:).^2);
end

w = vertcat(w{:});
g = vertcat(g{:}); 

costTorque = Function('costT', {w}, {eTorque});
costFiber = Function('costF', {w}, {eFiber});
costPennation = Function('costP', {w}, {ePennation});
cost = Function('cost', {w}, {J});




% % Create an NLP solver
%prob = struct('f', J, 'x', w1, 'g',g1);

% nlp prob : 
% "x" opt parameters, 'f' function to minimized, 'g' contraint function 
prob = struct('x', w, 'f', J , 'g',g); 

solver = nlpsol('solver', 'ipopt', prob);

 
% Solve the NLP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
    'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);

param_opt = w_opt(1:nparam);

