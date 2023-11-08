function [normalizedMuscleActiveForceLength] = MuscleForceLengthCSI(normalizedFiberLength,nMuscles)
% MUSCLEFORCELENGTHCSI Computes muscle force as a function of normalized fiber length
%
% Syntax:
%   muscleForce = MuscleForceLengthCSI(normalizedFiberLength, nMuscles)
%
% Description:
%   This function computes muscle force based on a curve created by cubic
%   spline interpolation of the points defined on the Gordon et al. curve,
%   as described in Buchanan's article.
%
% Input:
%   - normalizedFiberLength: SX variable of normalized fiber lengths.
%   - nMuscles: The number of muscles (type : double).
%
% Output:
%   - muscleForce: A SX function of muscle forces corresponding to the input
%     normalized fiber lengths.
% 
% Sources
% Gordon, A. M., Huxley, A. F., & Julian, F. J. (1966).
%     The variation in isometric tension with sarcomere length in 
%     vertebrate muscle fibres. The Journal of Physiology, 184(1), 170‑192.
%     https://doi.org/10.1113/jphysiol.1966.sp007909
%
% Buchanan, T. S., Lloyd, D. G., Manal, K., & Besier, T. F. (2004). 
%     Neuromusculoskeletal Modeling : Estimation of Muscle Forces and Joint 
%     Moments and Movements from Measurements of Neural Command. Journal of 
%     Applied Biomechanics, 20(4), 367‑395. 
%     https://doi.org/10.1123/jab.20.4.367
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data from Gordon and al. (1966).
optimalSarcomereLength = 2.05 ; 
sarcomereLength  =  [1.27, 1.67, 2.0, 2.25, 3.65] ;
y =  [0, 0.84, 1, 1, 0] ; 
x = sarcomereLength./ optimalSarcomereLength ; 

% Muscles can be assumed to produce zero force when shorter than 50% and 
% longer than 150% of their optimal length, as shown here. 
% Buchanan et al.(2004). 
x(end) = 1.5 ; 

% % figure of the cubic spline 
% figure("Name","Sarcomere active force length relationship (Normalized")
% xx = linspace(0,2,100);
% yy = spline(x,y,xx); % cubic spine 
% plot(x,y,'ok',xx,yy,'r')
% ylim([0 1.5])
% xlim([0 1.5])
% legend('data point','cubic spline')

spl = spline(x,y) ; % cubic spine 

% normalizedFiberLength = SX.sym('Fiber_length', nMuscles);
expression = normalizedFiberLength .* 0 ; 

for i = 1 : spl.pieces
    % Extract coefficients for the i-th pieces
    a = spl.coefs(i,4);
    b = spl.coefs(i,3);
    c = spl.coefs(i,2);
    d = spl.coefs(i,1);
    
    % Define the cubic polynomial for the i-th pieces
    poly_i = a + b*(normalizedFiberLength - x(i)) + c*(normalizedFiberLength - x(i)).^2 + d*(normalizedFiberLength - x(i)).^3 ;
 
    % Add the polynomial for the i-th segment to the expression using if_else
     segment_expression = if_else(normalizedFiberLength >= x(i), if_else(normalizedFiberLength < x(i+1), poly_i, 0), 0);

     expression = expression + segment_expression;
end

normalizedMuscleActiveForceLength = expression ; 

%% test of the function 
% % You may import casadi
% % import casadi 
% try 
% addpath('/Users/mickaelbegon/Downloads/casadi-3.6.3-osx64-matlab2018b/')
% catch 
% 
% end
% try
%  addpath('C:\Users\Stage\Desktop\Neuromusculoskeletal Modeling\Casadi')
% catch
% 
% end 
% 
% import casadi.*
% Create a numerical function for evaluation
% eval_func = Function('eval_func', {normalizedFiberLength(1)}, {normalizedMuscleActiveForceLength(1)});
% 
% % Generate x values in the desired range
% x_values = linspace(0, 1.6, 100);
% 
% % Evaluate the function at the specified x values
% y_values = zeros(size(x_values));
% for i = 1:length(x_values)
%     y_values(i) = full(eval_func(x_values(i)));
% end
% 
% % Plot the result
% figure("Name","Sarcomere active force length relationship (Normalized")
% plot(x_values, y_values);
% xlabel('normalizedFiberLength');
% ylabel('Muscle Passive Force');
% xlim([0, 1.6]);
% xlabel('Fiber length normalized ')
% ylabel('Fiber force normalized ')
end


