function [param_opt] = NLP_identification(known_parameters_num,muscle_tendon_parameters_num,unknown_parameters,casadiFun,Data,opts,ErrorIntialGuess)
% add information about the NLP 
    %% 1. Set Up function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% import Casadi
[~, name] = system('hostname');
if strcmp(name(1:9), '942-27984')
    addpath('C:\Users\Stage\Desktop\Doctorat\Manip_Neuromusculoskeletal_Modeling\Casadi')
elseif strcmp(name(1:8), '151302-1')
    addpath('C:\Users\amariani\Desktop\Thèse\Manip_Neuromusculoskeletal_Modeling\Casadi')
elseif strcmp(name(1:27), 'MacBook-Air-de-mickaelbegon')
    addpath('/Users/mickaelbegon/Downloads/casadi-3.6.3-osx64-matlab2018b/')
end
import casadi.*



nMuscles = 3 ; % Number of muscle in our model [TibialisAnterior, Soleus, Gastrocnemius]

rng(2,"twister") ; % 'seed' for randomisation

ntrials = 40 ; % Number of trials selected for the estimation of muscle tendon parameters
% ntrials = 18 ; % Number of trials selected for the estimation of muscle tendon parameters

    %% 2. Selection of trials 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2.1 Interst
opts = string(opts) ; 
if opts == "RANDOM"  % (Random seletion)
    selection = [1 randperm(size(Data,1), ntrials)];
    selection = sort(unique(selection));
    Data = Data(selection,:); % selected Data
    ntrials = size(Data,1);
elseif opts == "CHOSEN"
    ntrials = size(Data,1);
else
    error('Error in option selection : input option must be "RANDOM" or "CHOSEN" --> RANDOM is selected')
end

rng(2,"twister") ; % 'seed' for randomisation


    %% 3. NLP conficuration 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3.1 NLP Set Up 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3.1.1 Start with an empty NLP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

eTorque = [];
eFiber = [];
ePennation = [];

% 3.1.2 create muscle parameters variables (we should find muscle_tendon_parameters_num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = { w{:}, UP};
w0 =  [w0; muscle_tendon_parameters_num' .* ErrorIntialGuess' ]; % intial guess
lbw = [lbw; muscle_tendon_parameters_num' * 0.3]; % lower bound of variable
ubw = [ubw; muscle_tendon_parameters_num' * 3]; % upper bound of variable

% 3.1.3 Weightings in cost function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W_torque = 1; % 1 Nm
W_length = 0.005; % 1Nm correspond à 5 mm
W_angle = ((3/180) * pi); % 1 Nm correspond à 3 deg

%% 3.2 NLP equation creation (objective and cost function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for trial = 1:ntrials % for 1 to nb trials

    % extract experimental input data (mesured data)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data = Data(trial,:); % known variables
    a_trial = data(4:6) ; % muscle activation during the trial
    q_trial = [0, 0, 0, 0, data(2:3)] ; % skeleton configuration during the trial
    musculoskeletal_states_trial = [q_trial, known_parameters_num ] ; % muscleskeleton configuration during the trial

    % ℓmt : Muscle tendon unit length 
    UMTlength = full(casadiFun.getUMTLength(musculoskeletal_states_trial)) ;

    % define decision variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    str_trial = num2str(trial);
    fiberLength_k = SX.sym(['Fiber_length_' str_trial], nMuscles);
    tendonLengthening_k = SX.sym(['Tendon_Lengthening_' str_trial], nMuscles);
    PennationAngle_k = SX.sym(['Pennation_Angle_' str_trial], nMuscles);
    FT_k = SX.sym(['Tendon_Force_' str_trial], nMuscles);
    FM_k = SX.sym(['Muscle_Force_' str_trial], nMuscles);

    w0_k = data([19:21, 22:24, 25:27, 10:12, 13:15]) ; % mesured variables
    % w0_k = w0_k * ErrorInMesure(trial) ; % add noise in mesured variables
    w_k =  vertcat(FT_k,FM_k,tendonLengthening_k,fiberLength_k,PennationAngle_k);

    w = { w{:}, w_k}; % better to use tendon length
    w0 =  [w0; w0_k'];
    lbw = [lbw; w0_k' *.1]; % lower bound of variable
    ubw = [ubw; w0_k' * 3];

    K = vertcat(a_trial', UMTlength, UP) ;

    % Muscle tendon equilibium function  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    constraints = casadiFun.equilibriumError(w_k,K) ; % constraints function

    g = { g{:}, constraints};
    lbg = [lbg; zeros(15, 1)]; % lower bound of constraints
    ubg = [ubg; zeros(15, 1)]; % upper bound of constraints

    casadiFun.equilibriumError(w0_k, [a_trial'; UMTlength; muscle_tendon_parameters_num']); % equilibrium
    Torque_simulated = casadiFun.getJointMoment2( musculoskeletal_states_trial, FT_k) ; % simulated torque 

    % objectives functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    e_Torque = data(1) - Torque_simulated;
    e_Fiber = data(10:12)' - fiberLength_k;
    e_Pennation = data(13:15)' - PennationAngle_k;

    J = J + W_torque * e_Torque.^2; %add error on joint torque
    J = J + W_length * sum(e_Fiber.^2);% add error on tendon length and pennation angle
    J = J + W_angle * sum(e_Pennation.^2);

    eTorque = [eTorque; e_Torque]; % err between mesured and estimated torque (cost function)
    eFiber = [eFiber; e_Fiber]; % err between mesured and estimated fiber length (cost function)
    ePennation = [ePennation; e_Pennation]; % err between mesured and estimated pennation angle (cost function)
end

w = vertcat(w{:});
g = vertcat(g{:});


%% 3.3 NLP Solver 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "x" opt parameters, 'f' function to minimized, 'g' contraint function
prob = struct('x', vertcat(w), 'f', J, 'g', g); % problem to solve 
solver = nlpsol('solver', 'ipopt', prob); % solver set up --> type (ipopt) 

% Solve the NLP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
    'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x); 

param_opt = w_opt(1:length(UP)); % optimized Muscle Tendon Parameters (ℓom, φo, Fom, ℓst)  

% display between Reel and Estimated
disp('Error between Reel and Estimated Muscle-Tendon Parameters :')
disp(['Number of trials : ', num2str(ntrials)])
disp(['err ℓom : ', num2str(param_opt(1:3)')])
disp(['err φo : ', num2str(param_opt(4:6)')])
disp(['err Fom : ', num2str(param_opt(7:9)')])
disp(['err ℓst  : ', num2str(param_opt(10:12)')])
end