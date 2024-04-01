function plot_Muscle(representationFun)
%% set map 
x = 0.3 : 0.01 : 1.5 ; % fiber length
y = 0 : 0.05 : 1 ; % muscle activation
[X,Y] = meshgrid(x,y);
% Evaluate the equation at each point on the grid
Z = full(representationFun.normalized.fiber_total(Y,X));
Z2 =  full(representationFun.normalized.fiber_active(Y,X));
Z3 = full(representationFun.normalized.fiber_passive(X));

%% fig 1
figure ("Name","Fiber Force-Length-Activation relationship component","Color",[1 1 1])
surf(X, Y, Z2, 'FaceColor', 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.9);
hold on 
surf(X, Y, Z3, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.9);

% Add labels and title
xlabel('Fiber noramlized length','FontWeight','bold');
ylabel('Fiber normalized activation','FontWeight','bold');
zlabel('Fiber noramlized force','FontWeight','bold');
title('Fiber Force-Length-Activation relationship component');
xlim([0.3 1.5])
ylim([0 1])
legend('active fiber force','passive fiber force')


%% fig 2
% Plot the surface
figure("Name","Fiber Force-Length-Activation relationship","Color",[1 1 1])
%surf(X, Y, Z, 'FaceColor', 'interp', 'EdgeColor', 'interp', 'FaceAlpha', 0.9);
surf(X, Y, Z, 'FaceColor', 'interp', 'FaceAlpha', 0.9);

% Add labels and title
xlabel('Fiber noramlized length','FontWeight','bold');
ylabel('Fiber normalized activation','FontWeight','bold');
zlabel('Fiber noramlized force','FontWeight','bold');
title('Fiber Force-Length-Activation relationship');

xlim([0.3 1.5])
ylim([0 1])
%legend('Total fiber force')

end