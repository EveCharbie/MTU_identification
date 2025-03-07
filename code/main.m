                        %% Main 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up 
clear
close all
clc

    %% Load Neuromusculoskeletal :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import from Opensim Musculoskeletal Geometry (segment length, muscle insertion and via point) and Generic Muscle Tendon Parameters (ℓom, φo, Fom, ℓst)
[skeleton_num,muscle_tendon_parameters_num] = Model_OSIM2Mat() ; 

% Import Muscle Contraction Dynamics (muscle tendon equation from De
% Groote) --> Hill type model --> Fmt = f(a, ℓmt, νmt; Fom, ℓom, ℓst, φo).
% Note that in our model we ignore :
%       - fiber contraction velocity (νmt = 1)
%       - and electromechanical delay (a(t) = e(t)) 
[casadiFun,vizualizationFun,unknown_parameters,definition] = DeGrooteFunction() ; 

% test the current model (vizualization and how to use function 
testModel(skeleton_num,muscle_tendon_parameters_num,casadiFun,vizualizationFun)

    %% Data generator 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generation of Hypothetical Data with the model 
[Data,header] = HypotheticalDataGenerator(skeleton_num,muscle_tendon_parameters_num,casadiFun) ; 

    %% NLP  NonLinear Programming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimisation of Generic Muscle Tendon Parameters with experimental data 
load('Data.mat')
%% 1. random data ( to do : our protocol)
errParameters = NLP_identification_simulation(skeleton_num,muscle_tendon_parameters_num,unknown_parameters,casadiFun,Data,"RANDOM");

%% 2. our protocol
[ourData] = ourTrials(Data);

errParameters = NLP_identification_simulation(skeleton_num,muscle_tendon_parameters_num,unknown_parameters,casadiFun,ourData,"CHOSEN");


%% 2. test on real data (multistart)
% selection des données (protocol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sols = [];
% multistart
for i = 1 : 10
    % error in parameter estimation
    errEstimationParameter = 20/200 ; % error in percent (ex : 50/100 = 50% error)
    random_values = (errEstimationParameter/2) * randn(1, size(muscle_tendon_parameters_num,2)) + 1; % Generate random values from a normal distribution
    random_values(random_values < 1-errEstimationParameter) = 1-errEstimationParameter;
    random_values(random_values > 1+errEstimationParameter) = 1+errEstimationParameter; % between 50 % and 150%
    ErrorIntialGuess = random_values ;

    [param_opt] = NLP_identification(skeleton_num,muscle_tendon_parameters_num,unknown_parameters,casadiFun,NewData,"CHOSEN",ErrorIntialGuess);

    sols(i,:) = param_opt';
end