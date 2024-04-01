function [representationFun] = RepresentationFunction()
% import Casadi 
[~, name] = system('hostname');
if strcmp(name(1:9), '942-27984')
    addpath('C:\Users\Stage\Desktop\Doctorat\Manip_Neuromusculoskeletal_Modeling\Casadi')
elseif strcmp(name(1:27), 'MacBook-Air-de-mickaelbegon')
    addpath('/Users/mickaelbegon/Downloads/casadi-3.6.3-osx64-matlab2018b/')
end
import casadi.*

%% set up
a = SX.sym('muscle_activation');
lf = SX.sym('fiber_length');
lt = SX.sym('tendon_length');

%% fiber function only noramlized 
% Active force force-length (S2) 
b11 = 0.814483478343008 ; b21 = 1.055033428970575 ; b31 = 0.162384573599574 ; b41 = 0.063303448465465 ; % first Gaussian coefficents
b12 = 0.433004984392647 ; b22 = 0.716775413397760; b32 = -0.029947116970696 ; b42 = 0.200356847296188 ; % second Gaussian coefficents
b13 = 0.100 ; b23 = 1.000 ; b33 = 0.5 * sqrt(0.5) ; b43 = 0.000 ; % third Gaussian coefficents

num3 = lf - b23;
den3 = b33 + b43 * lf;
FMtilde3 = b13 * exp(-0.5 * num3.^2 ./ den3.^2);
num1 = lf - b21;
den1 = b31 + b41 * lf;
FMtilde1 = b11 * exp(-0.5 * num1.^2 ./ den1.^2);
num2 = lf - b22;
den2 = b32 + b42 * lf;
FMtilde2 = b12 * exp(-0.5 * num2.^2 ./ den2.^2);

Fiber_Active_F_L = a .* (FMtilde1 + FMtilde2 + FMtilde3); % Normalized equation 

% Passive fiber force-length (S3)
kpe = 4.0 ; e0 = 0.6 ; % fiber parameters

normalizedMusclePassiveForcePart1 =  0 ;
normalizedMusclePassiveForcePart2 = (exp(((kpe .* (lf - 1))./e0)) - 1)./ (exp(kpe) - 1) ;
Fiber_Passive_F_L = if_else(lf < 1, ...
    normalizedMusclePassiveForcePart1, ...
    normalizedMusclePassiveForcePart2); % if normalized length under 0 the force = 0            % Normalized equation 

% Total fiber force-length
Fiber_F_L = Fiber_Active_F_L + Fiber_Passive_F_L;


%% tendon 
% Tendon force-length (S1)
kT = 35; c1 = 0.200; c2 = 0.995; c3 = 0.250; % tendon parameters

Tendon_F_L = c1 .* exp(kT .* (lt - c2)) - c3; % Normalized equation 

%% function 
fiber_passive = Function('fiber_passive', ...
    {lf}, {Fiber_Passive_F_L}, ...
    {'noramlized fiber length'}, {'Fiber_Passive_Force'});

fiber_active = Function('fiber_active', ...
    {a,lf}, {Fiber_Active_F_L}, ...
    {'muscle_activity','noramlized fiber length'}, {'Fiber_Active_Force'});

fiber_total = Function('fiber_total', ...
    {a,lf}, {Fiber_F_L}, ...
    {'muscle_activity','noramlized fiber length'}, {'Fiber_Total_Force'});

tendon = Function('tendon', ...
    {lt}, {Tendon_F_L}, ...
    {'noramlized tendon length'}, {'Tendon_F_L'});

%% structure 
representationFun.normalized = struct(...
    'fiber_passive',fiber_passive,...
    'fiber_active', fiber_active,...
    'fiber_total', fiber_total,...
    'tendon', tendon) ; 

end 