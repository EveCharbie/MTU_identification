function [x0] = bestStart(a,lmtu,MTUparameters) 
% x0_1 : tendon force 
% x0_2 : muscle force 
% x0_3 : tendon lengthening  
% x0_4 : fiber length
% x0_5 : penation angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% muscle tendon parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%
Optimal_fiber_length = MTUparameters(1);
Pennation_angle_at_muscle_optimal_fiber_length = MTUparameters(2);
Maximal_isometric_muscle_force = MTUparameters(3);
Tendon_slack_length = MTUparameters(4);

% tendon lengthening
%%%%%%%%%%%%%%%%%%%%%
tendon_lengtheningRange = [0, .05];
x0_3 = randInRange(tendon_lengtheningRange,rand(1));
tendonLength = Tendon_slack_length + x0_3; 

if tendonLength<Tendon_slack_length
    tendonLength= Tendon_slack_length;
    x0_3 = 0;
end

% fiber length
%%%%%%%%%%%%%%%%%%%%%
fiberRange = [Optimal_fiber_length - Optimal_fiber_length*0.4,...
    Optimal_fiber_length + Optimal_fiber_length*0.4];            % en mètres
x0_4 = randInRange(fiberRange,rand(1));

% penation angle 
%%%%%%%%%%%%%%%%%%%%%
x0_5 =acos((lmtu - tendonLength) / x0_4);

% muscle force 
%%%%%%%%%%%%%%%%%%%%%
x0_2 = a * Maximal_isometric_muscle_force * (1 - abs(x0_4 - Optimal_fiber_length)/Optimal_fiber_length);

% tendon force 
%%%%%%%%%%%%%%%%%%%%%
x0_1 = a * Maximal_isometric_muscle_force * cos(x0_5);

% unknown  = vertcat(FT, FM, tendonLengthening, fiberLength, pennationAngle) ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = [x0_1,x0_2,x0_3,x0_4,x0_5];

if all(~isreal(x0))
    x0 = [Maximal_isometric_muscle_force * cos(Pennation_angle_at_muscle_optimal_fiber_length),...
        Maximal_isometric_muscle_force,...
        0.03,...
        Optimal_fiber_length,...
        Pennation_angle_at_muscle_optimal_fiber_length];
end
end 