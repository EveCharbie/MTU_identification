function plotMuscleForcelength(vizualizationFun,Optimal_fiber_length,Maximal_isometric_muscle_force,Muscle_activation,Fiber_length)
if nargin == 1 
    Optimal_fiber_length = 1; 
    Maximal_isometric_muscle_force = 1; 
    Muscle_activation = NaN(1);
    Fiber_length = NaN(1);
elseif nargin == 3
    Muscle_activation = NaN(1);
    Fiber_length = NaN(1);
end

%%
% about activation
activation = linspace(0,1,100); % 100 valeurs d'activation

% about fiber length
fiber_normalized_length = linspace(0.2,1.6,100); % 100 valeurs de longueur
fiber_length = fiber_normalized_length* Optimal_fiber_length;
fiber_current_normalized_length = Fiber_length/Optimal_fiber_length;

% about force 
[ACT, LENGTHNORMALIZED] = meshgrid(activation, fiber_normalized_length);
[~, LENGTH] = meshgrid(activation, fiber_length);


    % normalized
NormalizedActiveForce = ACT .* full(vizualizationFun.getNormalizedMuscleActiveForce(LENGTHNORMALIZED,1));
NormalizedPassiveForce = full(vizualizationFun.getNormalizedMusclePassiveForce(LENGTHNORMALIZED,1));
NormalizedForce = full(vizualizationFun.getNormalizedMuscleForce(ACT,LENGTHNORMALIZED,1));

NormalizedActiveForceCurrent = Muscle_activation .* full(vizualizationFun.getNormalizedMuscleActiveForce(fiber_current_normalized_length,1));
NormalizedPassiveForceCurrent = full(vizualizationFun.getNormalizedMusclePassiveForce(fiber_current_normalized_length,1));
NormalizedForceCurrent = full(vizualizationFun.getNormalizedMuscleForce(Muscle_activation,fiber_current_normalized_length,1));

    % no-normalized
activeForce = NormalizedActiveForce .* Maximal_isometric_muscle_force;
passiveForce = NormalizedPassiveForce .* Maximal_isometric_muscle_force;
Force = NormalizedForce .* Maximal_isometric_muscle_force;

activeForceCurrent = NormalizedActiveForceCurrent .* Maximal_isometric_muscle_force;
passiveForceCurrent = NormalizedPassiveForceCurrent .* Maximal_isometric_muscle_force;
ForceCurrent = NormalizedForceCurrent .* Maximal_isometric_muscle_force;

%% plot 
%% normalized
% Affichage du graphique 3D
figure("Name","Muscle force-length relationship")
subplot(3,2,1)
hold on 
plot3(fiber_current_normalized_length,Muscle_activation,NormalizedActiveForceCurrent,'+r')
surf(LENGTHNORMALIZED, ACT, NormalizedActiveForce);
colorbar; % Pour une meilleure visualisation
shading interp; % Lissage du rendu
xlabel('Normalized fiber length','FontWeight','bold')
ylabel('Normalized fiber activation','FontWeight','bold')
zlabel('Normalized fiber force','FontWeight','bold')
title('Normalized active force-lenght relationship');
xlim([min(LENGTHNORMALIZED,[],"all") max(LENGTHNORMALIZED,[],"all")])
ylim([min(ACT,[],"all") max(ACT,[],"all")])
zlim([min(NormalizedActiveForce,[],"all") max(NormalizedActiveForce,[],"all")])
view(-5,25)
hold off

subplot(3,2,3)
hold on 
plot3(fiber_current_normalized_length,Muscle_activation,NormalizedPassiveForceCurrent,'+r')
surf(LENGTHNORMALIZED, ACT, NormalizedPassiveForce);
colorbar; % Pour une meilleure visualisation
shading interp; % Lissage du rendu
xlabel('Normalized fiber length','FontWeight','bold')
ylabel('Normalized fiber activation','FontWeight','bold')
zlabel('Normalized fiber force','FontWeight','bold')
title('Normalized passive force-lenght relationship');
xlim([min(LENGTHNORMALIZED,[],"all") max(LENGTHNORMALIZED,[],"all")])
ylim([min(ACT,[],"all") max(ACT,[],"all")])
zlim([min(NormalizedPassiveForce,[],"all") max(NormalizedPassiveForce,[],"all")])
view(-5,25)
hold off

subplot(3,2,5)
hold on 
plot3(fiber_current_normalized_length,Muscle_activation,NormalizedForceCurrent,'+r')
surf(LENGTHNORMALIZED, ACT, NormalizedForce);
colorbar; % Pour une meilleure visualisation
shading interp; % Lissage du rendu
xlabel('Normalized fiber length','FontWeight','bold')
ylabel('Normalized fiber activation','FontWeight','bold')
zlabel('Normalized fiber force','FontWeight','bold')
title('Normalized total force-lenght relationship');
xlim([min(LENGTHNORMALIZED,[],"all") max(LENGTHNORMALIZED,[],"all")])
ylim([min(ACT,[],"all") max(ACT,[],"all")])
zlim([min(NormalizedForce,[],"all") max(NormalizedForce,[],"all")])
view(-5,25)
hold off



%% not normalized
% Affichage du graphique 3D
subplot(3,2,2)
hold on 
plot3(Fiber_length,Muscle_activation,activeForceCurrent,'+r') 
surf(LENGTH, ACT, activeForce);
colorbar; % Pour une meilleure visualisation
shading interp; % Lissage du rendu
xlabel('Fiber length','FontWeight','bold')
ylabel('Fiber activation','FontWeight','bold')
zlabel('Fiber force','FontWeight','bold')
title('Active force-lenght relationship');
xlim([min(LENGTH,[],"all") max(LENGTH,[],"all")])
ylim([min(ACT,[],"all") max(ACT,[],"all")])
zlim([min(activeForce,[],"all") max(activeForce,[],"all")])
view(-5,25)
hold off

subplot(3,2,4)
hold on 
plot3(Fiber_length,Muscle_activation,passiveForceCurrent,'+r')
surf(LENGTH, ACT, passiveForce);
colorbar; % Pour une meilleure visualisation
shading interp; % Lissage du rendu
xlabel('Fiber length','FontWeight','bold')
ylabel('Fiber activation','FontWeight','bold')
zlabel('Fiber force','FontWeight','bold')
title('Passive force-lenght relationship');
xlim([min(LENGTH,[],"all") max(LENGTH,[],"all")])
ylim([min(ACT,[],"all") max(ACT,[],"all")])
zlim([min(passiveForce,[],"all") max(passiveForce,[],"all")])
view(-5,25)
hold off

subplot(3,2,6)
hold on 
plot3(Fiber_length,Muscle_activation,ForceCurrent,'+r') 
surf(LENGTH, ACT, Force);
colorbar; % Pour une meilleure visualisation
shading interp; % Lissage du rendu
xlabel('Fiber length','FontWeight','bold')
ylabel('Fiber activation','FontWeight','bold')
zlabel('Fiber force','FontWeight','bold')
title('Total force-lenght relationship');
xlim([min(LENGTH,[],"all") max(LENGTH,[],"all")])
ylim([min(ACT,[],"all") max(ACT,[],"all")])
zlim([min(Force,[],"all") max(Force,[],"all")])
view(-5,25)
hold off

end 