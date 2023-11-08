function [] = Force_length_representation(all_states_num,muscle_tendon_parameters_num,casadiFun,NameMuscles)
%% errrror 
%% forces illustration function 
% it shows force length relation ship of all muscle and tendon 
nMuscles = length(NameMuscles) ; 
plotDefinition = 100 ; 
Length_Muscle_Relationship = linspace(0.2, 1.4, plotDefinition) ;  
Lengthing_Tendon_Relationship = linspace(-0.05, .05, plotDefinition) ; 
Length_Tendon_Relationship = Lengthing_Tendon_Relationship +1 ; 

l0m_num = muscle_tendon_parameters_num(1:3);
phi0_num = muscle_tendon_parameters_num(4:6);
f0m_num =  muscle_tendon_parameters_num(7:9);
lst_num =  muscle_tendon_parameters_num(10:12);

Noramized_Fiber_Length_num = full(casadiFun.normalizeFiberLength(all_states_num,muscle_tendon_parameters_num)) ; 
Noramized_Tendon_Length_num = full(casadiFun.normalizeTendonLength(all_states_num,muscle_tendon_parameters_num)) ; 

% initialisation (for perf)
Active_Muscle_Force_Relationship = nan(1,plotDefinition) ; 
Passive_Muscle_Force_Relationship = nan(1,plotDefinition) ; 
Muscle_Force_Relationship = nan(1,plotDefinition) ; 
Tendon_Force_Relationship = nan(1,plotDefinition) ; 

%% relationship muscle by muscle 
all_states_figure = full(all_states_num) ; 
for ii = 1 : plotDefinition
    all_states_figure(44: 46) = ones(1,3) .* Lengthing_Tendon_Relationship(ii) ;
    all_states_figure(41: 43) = ones(1,3) .* Length_Tendon_Relationship(ii) ;

    %temp  = full(casadiFun.getTendonForce(all_states_figure,muscle_tendon_parameters_num))
    Tendon_Force_Relationship(:,ii) = full(casadiFun.getTendonForce(all_states_figure,muscle_tendon_parameters_num))

    Active_Muscle_Force_Relationship(:,ii)  = full(casadiFun.getMuscleForce(all_states_figure,muscle_tendon_parameters_num)) ; 
    all_states_figure(end-nMuscles - nMuscles+1 : end-nMuscles) ;
end
% Generate x values in the desired range
x_values = linspace(-.05, .05, 100);



for ii = 1 : plotDefinition
    all_states_figure(end-nMuscles - nMuscles+1 : end-nMuscles) = [Lengthing_Tendon_Relationship(ii), Lengthing_Tendon_Relationship(ii),Lengthing_Tendon_Relationship(ii)] ; 


end
% Evaluate the function at the specified x values
y_values = zeros(size(x_values));
for i = 1:length(x_values)
    y_values(i) = full(eval_func(x_values(i),0.35));
end

% Plot the result
figure;
plot(x_values, y_values);
xlabel('normalizedFiberLength');
ylabel('Muscle Passive Force');
xlim([-.05 .05]);
grid on;

















for ii = 1 : plotDefinition
    all_states_figure(end-nMuscles - nMuscles+1 : end-nMuscles) = [Lengthing_Tendon_Relationship(ii), Lengthing_Tendon_Relationship(ii),Lengthing_Tendon_Relationship(ii)] ; 
    all_states_figure(end-(nMuscles*2) - nMuscles+1 : end-(nMuscles*2)) = [Length_Muscle_Relationship(ii), Length_Muscle_Relationship(ii),Length_Muscle_Relationship(ii)] ; 

    Active_Muscle_Force_Relationship(:,ii) = full(casadiFun.getMuscleActiveForce(all_states_figure,muscle_tendon_parameters_num)) ; 
    Passive_Muscle_Force_Relationship(:,ii) = full(casadiFun.getMusclePassiveForce(all_states_figure,muscle_tendon_parameters_num)) ; 
    Muscle_Force_Relationship(:,ii) = Active_Muscle_Force_Relationship(:,ii) + Passive_Muscle_Force_Relationship(:,ii);
    Tendon_Force_Relationship(:,ii) = full(casadiFun.getTendonForce(all_states_figure,muscle_tendon_parameters_num)) ; 
end


Active_Muscle_Force = full(casadiFun.getMuscleActiveForce(all_states_num,muscle_tendon_parameters_num)) ;
Passive_Muscle_Force= full(casadiFun.getMusclePassiveForce(all_states_num,muscle_tendon_parameters_num)) ;
Muscle_Force =  Active_Muscle_Force + Passive_Muscle_Force;
Tendon_Force = full(casadiFun.getTendonForce(all_states_num,muscle_tendon_parameters_num)) ;

%% plots of muscle forces 
figure ('Name','Mucles forces representation')
for i = 1 : nMuscles
subplot(1,nMuscles,i)
plot(Length_Muscle_Relationship, Active_Muscle_Force_Relationship(i,:) ,':k')
hold on 
plot(Length_Muscle_Relationship, Passive_Muscle_Force_Relationship(i,:) ,':k')
plot(Length_Muscle_Relationship, Muscle_Force_Relationship(i,:) ,'k')

plot(Noramized_Fiber_Length_num(i),Active_Muscle_Force(i),'or')
plot(Noramized_Fiber_Length_num(i),Passive_Muscle_Force(i),'or')
plot(Noramized_Fiber_Length_num(i),Muscle_Force(i),'or')
legend([{'Active muscle force [CC]'},{'Passive muscle force [CEP]'},{'Total muscle force [CC + CEP]'}])
%title([NameMuscles{i}, '- EMG : ', num2str(a_num(i))])
title(NameMuscles{i})

xlabel('Normalized lenght', 'FontWeight','bold')
ylabel('Force (N)', 'FontWeight','bold')
grid on 
end

% plot of tendon force 
figure ('Name','Tendons forces representation')
for i =  1 : nMuscles
subplot(1,nMuscles,i)
plot(Length_Tendon_Relationship, Tendon_Force_Relationship(i,:),':k')
hold on 
plot(Noramized_Tendon_Length_num(i), Tendon_Force(i),'or')
%title([NameMuscles{i}, '- EMG : ', num2str(a_num(i))])
title(NameMuscles{i})

xlabel('Normalized lenght', 'FontWeight','bold')
ylabel('Force (N)', 'FontWeight','bold')
xlim([0.95 1.05])
grid on 
end

PC_Max_Force = nan(nMuscles) ;
% muscle force in percent of max muscle force
for i = 1 : nMuscles
PC_Max_Force(i) = (Active_Muscle_Force(i) / f0m_num(i)) *100 ; 
end

% print in command window
fprintf('Muscle are resprectively at \n ')
for i = 1 : nMuscles
fprintf([ '   - ', NameMuscles{i}, ' : ' , num2str(PC_Max_Force(i)) , ' percent of his maximal force. \n'])
end
fprintf('\n')

fprintf('Muscle normalized length are resprectively at \n ')
for i = 1 : nMuscles
fprintf([ '   - ', NameMuscles{i}, ' : ' , num2str(Noramized_Fiber_Length_num(i)) , ' . \n'])
end
fprintf('\n')

fprintf('Tendon normalized length are resprectively at \n ')
for i = 1 : nMuscles
fprintf([ '   - ', NameMuscles{i}, ' : ' , num2str(Noramized_Tendon_Length_num(i)) , ' . \n'])
end

end
