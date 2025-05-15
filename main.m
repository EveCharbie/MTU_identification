                        %% Main 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up 
clear
close all
clc

current_folder = cd; 

folders = struct( ...
    'main', current_folder, ...
    'fun', fullfile(current_folder, 'matFun'), ...
    'osim_model', fullfile(current_folder, 'osimModel'), ...
    'simulatedData', fullfile(current_folder, 'simulatedData') ...
);
    %% Load Neuromusculoskeletal :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import from Opensim Musculoskeletal Geometry and Generic Muscle Tendon Parameters (ℓom, φo, Fom, ℓst)
[known_parameters_num,muscle_tendon_parameters_num] = Model_OSIM2Mat(folders.osim_model) ; 

addpath(folders.fun)

% Import Muscle Contraction Dynamics (muscle tendon equation from De
% Groote) --> Hill type model --> Fmt = f(a, ℓmt, νmt; Fom, ℓom, ℓst, φo).
% Note that in our model we ignore :
%       - fiber contraction velocity (νmt = 1)
%       - and electromechanical delay (a(t) = e(t)) 
[casadiFun,unknown_parameters,definition] = DeGrooteFunction() ; 

[casadiFun2,unknown_parameters2,definition2] = DeGrooteFunctionTest();
       %% test the neuromusculo model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
testModel(known_parameters_num,muscle_tendon_parameters_num,casadiFun)

    %% Data generator 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generation of Hypothetical Data with the model 
[Data,header] = HypotheticalDataGenerator(known_parameters_num,muscle_tendon_parameters_num,casadiFun) ; 

    %% NLP  NonLinear Programming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimisation of Generic Muscle Tendon Parameters with experimental data 
load(fullfile(folders.simulatedData,'Data.mat'))
%% 1. random data ( to do : our protocol)
errParameters = NLP_identification_simulation(known_parameters_num,muscle_tendon_parameters_num,unknown_parameters,casadiFun,Data,"RANDOM");

%% 2. test on real data (multistart) 
% selection des données (protocol) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Définition des valeurs à filtrer 
val_col3 = ([20, 10, 0, -10, -20] / 180) * pi; 
val_col4 = ([80,  0] / 180) * pi; 
val_col4 = val_col3
val_col456 = [0, 0.2, 0.4, 0.6]; 
 
% Création du masque de filtrage 
mask = (Data(:,2) == 0) & ... 
       ismember(Data(:,3), val_col3) & ... 
       (ismember(Data(:,4), val_col456) | ismember(Data(:,5), val_col456) | ismember(Data(:,6), val_col456)) & ... 
       ~( (Data(:,4) == 1) | (Data(:,5) == 1) | (Data(:,6) == 1) ); % Exclure les lignes contenant 1 
 
% mask2 = (Data(:,2) == max(Data(:,2))) & ... 
%        ismember(Data(:,3), val_col3([4,5])) & ... 
%        (ismember(Data(:,4), val_col456) | ismember(Data(:,5), val_col456) | ismember(Data(:,6), val_col456)) & ... 
%        ~( (Data(:,4) == 1) | (Data(:,5) == 1) | (Data(:,6) == 1) ); % Exclure les lignes contenant 1
mask2 = ismember(Data(:,2), val_col4) & ... 
       ismember(Data(:,3), val_col3) & ... 
       (ismember(Data(:,4), val_col456) | ismember(Data(:,5), val_col456) | ismember(Data(:,6), val_col456)) & ... 
       ~( (Data(:,4) == 1) | (Data(:,5) == 1) | (Data(:,6) == 1) ); % Exclure les lignes contenant 1

mask3 = mask+mask2 >= 1;

% Extraction des lignes correspondant aux conditions 
NewData = Data(mask3, :); 
sols = []; 
cost = [];

Noise = [0.2, 0.2, 0.40, 0.2]; % 
% multistart 
% rng('shuffle');
nTrials = 1000;
for i = 1 : nTrials
    rng(i,"twister") ; % 'seed' for randomisation
    intialGuess = randParameter(muscle_tendon_parameters_num,Noise) ; 
 
    [param_opt,cost(i)] = NLP_identification(known_parameters_num,muscle_tendon_parameters_num,unknown_parameters,casadiFun,NewData,"CHOSEN",intialGuess); 
 
    pStart(i,:) = intialGuess;
    sols(i,:) = param_opt'; 

    disp([num2str(i), ' out of ',num2str(nTrials)])
end


err_ = sols - muscle_tendon_parameters_num;
err_normalized = (err_ ./muscle_tendon_parameters_num)* 100;
cost = log(cost);


%% Results 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
figure('Name','results of musltistart','WindowState','maximized')
% tibialis 
subplot(3,4,1)
title('tibailis')
hold on
plot(abs(sols(:,1)),cost(:),'.k')
plot(abs(sols((cost<0.1),1)),cost(cost<0.1),'or')
xlabel('l0m', 'FontWeight','bold')
ylabel('log_C_o_s_t', 'FontWeight','bold')
hold off

subplot(3,4,2)
hold on
plot(abs(sols(:,4)),cost(:),'.k')
plot(abs(sols((cost<0.1),4)),cost(cost<0.1),'or')
xlabel('phi0', 'FontWeight','bold')
% ylabel('log_C_o_s_t', 'FontWeight','bold')
hold off

subplot(3,4,3)
hold on
plot(abs(sols(:,7)),cost(:),'.k')
plot(abs(sols((cost<0.1),7)),cost(cost<0.1),'or')
xlabel('f0m', 'FontWeight','bold')
% ylabel('log_C_o_s_t', 'FontWeight','bold')
hold off

subplot(3,4,4)
hold on
plot(abs(sols(:,10)),cost(:),'.k')
plot(abs(sols((cost<0.1),10)),cost(cost<0.1),'or')
xlabel('tsl', 'FontWeight','bold')
% ylabel('log_C_o_s_t', 'FontWeight','bold')
hold off


% soleus 
subplot(3,4,5)
title('soleus')
hold on
plot(abs(sols(:,2)),cost(:),'.k')
plot(abs(sols((cost<0.1),2)),cost(cost<0.1),'or')
xlabel('l0m', 'FontWeight','bold')
ylabel('log_C_o_s_t', 'FontWeight','bold')
hold off

subplot(3,4,6)
hold on
plot(abs(sols(:,5)),cost(:),'.k')
plot(abs(sols((cost<0.1),5)),cost(cost<0.1),'or')
xlabel('phi0', 'FontWeight','bold')
% ylabel('log_C_o_s_t', 'FontWeight','bold')
hold off

subplot(3,4,7)
hold on
plot(abs(sols(:,8)),cost(:),'.k')
plot(abs(sols((cost<0.1),8)),cost(cost<0.1),'or')
xlabel('f0m', 'FontWeight','bold')
% ylabel('log_C_o_s_t', 'FontWeight','bold')
hold off

subplot(3,4,8)
hold on
plot(abs(sols(:,11)),cost(:),'.k')
plot(abs(sols((cost<0.1),11)),cost(cost<0.1),'or')
xlabel('tsl', 'FontWeight','bold')
% ylabel('log_C_o_s_t', 'FontWeight','bold')
hold off

% gast 
subplot(3,4,9)
title('gastrocnemisue ')
hold on
plot(abs(sols(:,3)),cost(:),'.k')
plot(abs(sols((cost<0.1),3)),cost(cost<0.1),'or')
xlabel('l0m', 'FontWeight','bold')
ylabel('log_C_o_s_t', 'FontWeight','bold')
hold off

subplot(3,4,10)
hold on
plot(abs(sols(:,6)),cost(:),'.k')
plot(abs(sols((cost<0.1),6)),cost(cost<0.1),'or')
xlabel('phi0', 'FontWeight','bold')
% ylabel('log_C_o_s_t', 'FontWeight','bold')
hold off

subplot(3,4,11)
hold on
plot(abs(sols(:,9)),cost(:),'.k')
plot(abs(sols((cost<0.1),9)),cost(cost<0.1),'or')
xlabel('f0m', 'FontWeight','bold')
% ylabel('log_C_o_s_t', 'FontWeight','bold')
hold off

subplot(3,4,12)
hold on
plot(abs(sols(:,12)),cost(:),'.k')
plot(abs(sols((cost<0.1),12)),cost(cost<0.1),'or')
xlabel('tsl', 'FontWeight','bold')
% ylabel('log_C_o_s_t', 'FontWeight','bold')
hold off