function plot_Muscle_current(representationFun,muscleactivation,normalized_fiber_lenght,maximalIsometricForce)
%% set map 
x = 0.3 : 0.01 : 1.5 ; % fiber length
y = 0 : 0.05 : 1 ; % muscle activation
[X,Y] = meshgrid(x,y);
% Evaluate the equation at each point on the grid
Z = full(representationFun.normalized.fiber_total(Y,X)) * maximalIsometricForce;

xtrial = normalized_fiber_lenght; 
ytrial = muscleactivation ; 

ztrial = full(representationFun.normalized.fiber_total(ytrial,xtrial)) * maximalIsometricForce;
%% fig 2
% Plot the surface
fcolor = [0.6 0.6 0.6];
figure("Name","Fiber Force-Length-Activation relationship","Color",[1 1 1])
%surf(X, Y, Z, 'FaceColor', 'interp', 'EdgeColor', 'interp', 'FaceAlpha', 0.9);
surf(X, Y, Z, 'FaceColor', fcolor, 'FaceAlpha', 0.6);
hold on 

plot3(xtrial,ytrial,ztrial,'dr','MarkerSize',6,'LineWidth',4)

% Add labels and title
xlabel('Fiber noramlized length','FontWeight','bold');
ylabel('Fiber normalized activation','FontWeight','bold');
zlabel('Fiber noramlized force','FontWeight','bold');
title('Fiber Force-Length-Activation relationship');

xlim([0.3 1.5])
ylim([0 1])
%legend('Total fiber force')

end