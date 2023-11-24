clear, close all,

% import casadi
[ret, name] = system('hostname');
if strcmp(name(1:9), '942-27984')
    addpath('C:\Users\Stage\Desktop\Neuromusculoskeletal Modeling\Casadi')
elseif strcmp(name(1:27), 'MacBook-Air-de-mickaelbegon')
    addpath('/Users/mickaelbegon/Downloads/casadi-3.6.3-osx64-matlab2018b/')
end

import casadi.*

% equation of our model 
De_Groote_opensim ;

% musculoskeletal scale 
[known_parameters_num,muscle_tendon_parameters_num] = Opensim_extraction() ;



% Call data_generator
% hypothetical_data_generator ;
load('Data.mat')
ntrials = 50; %size(Data,1) ;
selection = [1 randperm(size(Data,1), ntrials)];
selection = sort(unique(selection));
Data = Data(selection,:);
ntrials = size(Data,1);


% load best start value 
load('StartEquilibrium.mat')
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

% % create muscle parameters variables
w = { w{:}, UP};
% % we should find muscle_tendon_parameters_num
w0 =  [w0; muscle_tendon_parameters_num' ]; %TODO: add noise
lbw = [lbw; muscle_tendon_parameters_num' * 0.5]; % lower bound of variable
ubw = [ubw; muscle_tendon_parameters_num' * 2]; % upper bound of variable


% header = header = {'Torque','q1','q2',...
%     'activation_tibialis','activation_soleus','activation_gastrocnemius',...
%     'fiber_tibialis','fiber_soleus','fiber_gastrocnemius',...
%     'phi_tibialis','phi_soleus','phi_gastrocnemius'} ;

%randomValues = randperm(maxValue, numValues);


% create variable for each trial

% Weightings in cost function
W_torque = 1; % 1 Nm
W_length = 1;% / 0.005; %todo 1Nm correspond à 5 mm
W_angle = 1;% / ((3/180) * pi);  % 1 Nm correspond à 3 deg
W_tendon = 1;% / 0.005; %todo 1Nm correspond à 5 mm


eTorque = [];
eFiber = [];
ePennation = [];
eTendon = [] ;

for trial = 1:ntrials % for 1 to nb trials

    % extract experimental input data
    data = Data(trial,:); % known variables 
    a_trial = data(4:6) ; % muscle activation during the trial
    q_trial = [0, 0, 0, 0, data(2:3)] ; % skeleton configuration during the trial
    musculoskeletal_states_trial = [q_trial, known_parameters_num ] ; % muscleskeleton configuration during the trial
    neuromusculoskeletal_states_trial = [a_trial, musculoskeletal_states_trial] ; % Neuromusculoskeletal states
    p_trial = vertcat(neuromusculoskeletal_states_trial', UP) ;

        % input 
    % muscle activation
    a_ta= a_trial(1);
    a_sol = a_trial(2);
    a_gast= a_trial(3);

    % UMT length
    UMTlength = full(casadiFun.getUMTLength(musculoskeletal_states_trial)) ;
    length_UMT_ta = UMTlength(1) ;
    length_UMT_sol = UMTlength(2) ;
    length_UMT_gast = UMTlength(3) ;
    

    % define decision variables
    str_trial = num2str(trial);
    fiberLength_k = SX.sym(['Fiber_length_' str_trial], nMuscles);
    tendonLengthening_k = SX.sym(['Tendon_Lengthening_' str_trial], nMuscles);
    PennationAngle_k = SX.sym(['Pennation_Angle_' str_trial], nMuscles);
    FT_k = SX.sym(['Tendon_Force_' str_trial], nMuscles);
    FM_k = SX.sym(['Muscle_Force_' str_trial], nMuscles);
        

    w0_k = data([19:21, 22:24, 25:27, 10:12, 13:15]);
    w_k =  vertcat(FT_k,FM_k,tendonLengthening_k,fiberLength_k,PennationAngle_k);

    w = { w{:}, w_k}; % better to use tendon length 
    w0 =  [w0; w0_k']; 
    lbw = [lbw; w0_k' *.1]; % lower bound of variable
    ubw = [ubw; w0_k' * 3]; 
    
    %known = vertcat(a, LUMT, muscleTendonParameters) ; 
    K = vertcat(a_trial', UMTlength, UP) ;

    % constraints
    constraints = casadiFun.equilibriumError1(w_k,K) ; 

    casadiFun.equilibriumError1(w0_k, [a_trial'; UMTlength; muscle_tendon_parameters_num'])
    
    
    g = { g{:}, constraints};
    lbg = [lbg; zeros(15, 1)]; % lower bound of constraints
    ubg = [ubg; zeros(15, 1)]; % upper bound of constraints

    Torque_simulated = casadiFun.getJointMoment2( musculoskeletal_states_trial, FT_k) ;
    
    % objective
    e_Torque = data(1) - Torque_simulated;
    e_Fiber = data(10:12)' - fiberLength_k;
    e_Pennation = data(13:15)' - PennationAngle_k;

    J = J + W_torque * e_Torque.^2; %add error on joint torque
    J = J + W_length * sum(e_Fiber.^2);% add error on tendon length and pennation angle
    J = J + W_angle * sum(e_Pennation.^2);

    eTorque = [eTorque; e_Torque];
    eFiber = [eFiber; e_Fiber];
    ePennation = [ePennation; e_Pennation];

end

w = vertcat(w{:});
g = vertcat(g{:}); 


costTorque = Function('costT', {w}, {eTorque});
costFiber = Function('costF', {w}, {eFiber});
costPennation = Function('costP', {w}, {ePennation});
cost = Function('cost', {w}, {J});
const = Function('const', {w}, {g});



% eval cost at w0
costTorque(w0)
costFiber(w0)
costPennation(w0)
const(w0)


% % Create an NLP solver
%prob = struct('f', J, 'x', w1, 'g',g1);

% nlp prob : 
% "x" opt parameters, 'f' function to minimized, 'g' contraint function 
prob = struct('x', vertcat(w), 'f', J, 'g', g); 

solver = nlpsol('solver', 'ipopt', prob);

% Solve the NLP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
    'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);

param_opt = w_opt(1:length(UP));


err = param_opt' - muscle_tendon_parameters_num ; 


