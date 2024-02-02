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
errParameters = NLP_identification_simulation(known_parameters_num,muscle_tendon_parameters_num,unknown_parameters,casadiFun,Data);