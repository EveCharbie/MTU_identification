                        %% Main 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up 
clear
close all
clc

    %% Load Neuromusculoskeletal :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import from Opensim Musculoskeletal Geometry and Generic Muscle Tendon Parameters (ℓom, φo, Fom, ℓst)
[known_parameters_num,muscle_tendon_parameters_num] = Model_OSIM2Mat() ; 

% Import Muscle Contraction Dynamics (muscle tendon equation from De
% Groote) --> Hill type model --> Fmt = f(a, ℓmt, νmt; Fom, ℓom, ℓst, φo).
% Note that in our model we ignore :
%       - fiber contraction velocity (νmt = 1)
%       - and electromechanical delay (a(t) = e(t)) 
[casadiFun,unknown_parameters,definition] = DeGrooteFunction() ; 

    %% Data generator 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generation of Hypothetical Data with the model 
[Data,header] = HypotheticalDataGenerator(known_parameters_num,muscle_tendon_parameters_num,casadiFun) ; 

    %% NLP  NonLinear Programming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimisation of Generic Muscle Tendon Parameters with experimental data 
load('Data.mat')
%% 1. random data ( to do : our protocol)
errParameters = NLP_identification_simulation(known_parameters_num,muscle_tendon_parameters_num,unknown_parameters,casadiFun,Data,"RANDOM");

%% 2. test on real data (multistart)
% selection des données (protocol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Définition des valeurs à filtrer
val_col3 = ([20, 10, 0, -10, -20] / 180) * pi;
val_col456 = [0, 0.2, 0.4, 0.6];

% Création du masque de filtrage
mask = (Data(:,2) == 0) & ...
       ismember(Data(:,3), val_col3) & ...
       (ismember(Data(:,4), val_col456) | ismember(Data(:,5), val_col456) | ismember(Data(:,6), val_col456)) & ...
       ~( (Data(:,4) == 1) | (Data(:,5) == 1) | (Data(:,6) == 1) ); % Exclure les lignes contenant 1

% Extraction des lignes correspondant aux conditions
NewData = Data(mask, :);
sols = [];
% multistart
for i = 1 : 10
    % error in parameter estimation
    errEstimationParameter = 50/100 ; % error in percent (ex : 50/100 = 50% error)
    random_values = (errEstimationParameter/2) * randn(1, size(muscle_tendon_parameters_num,2)) + 1; % Generate random values from a normal distribution
    random_values(random_values < 1-errEstimationParameter) = 1-errEstimationParameter;
    random_values(random_values > 1+errEstimationParameter) = 1+errEstimationParameter; % between 50 % and 150%
    ErrorIntialGuess = random_values ;

    [param_opt] = NLP_identification(known_parameters_num,muscle_tendon_parameters_num,unknown_parameters,casadiFun,NewData,"CHOSEN",ErrorIntialGuess)

    sols(i,:) = param_opt';
end