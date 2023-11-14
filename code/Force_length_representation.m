function [] = Force_length_representation(all_states_num,muscle_tendon_parameters_num,casadiFun,NameMuscles)
%% set up
nMuscles = length (NameMuscles) ; % nb muscle
definition = 200; %nb pts in the representation plot

% muscle parameters
l0m = muscle_tendon_parameters_num(1:3);
phi0 = muscle_tendon_parameters_num(4:6);
f0m =  muscle_tendon_parameters_num(7:9);
lst =  muscle_tendon_parameters_num(10:12);
a = full(all_states_num(end-2 : end)) ;

% structures length
fiberLength = full(casadiFun.normalizeFiberLength(all_states_num,muscle_tendon_parameters_num)) ;
tendonLength = full(casadiFun.normalizeTendonLength(all_states_num,muscle_tendon_parameters_num)) ;


% vecteur force length hypothetique
fiberLengthRepresentation = linspace(0, 1.6, definition);
tendonLengtheningRepresentation = linspace(-.05,.05, definition);
%casadiFun.

% Muscle
% passive force
for i = 1 : nMuscles
    for ii= i : definition
        MusclePassiveForceRepresentation(i,ii) = full(casadiFun.representationMusclePassiveForce(fiberLengthRepresentation(ii), 1 , f0m(i))) ;
    end
end
% active force
for i = 1 : nMuscles
    for ii= i : definition
        MuscleActiveForceRepresentation(i,ii) = full(casadiFun.representationMuscleActiveForceLength(a(i),fiberLengthRepresentation(ii), 1 , f0m(i))) ;
    end
end

MuscleForceRepresentation = MusclePassiveForceRepresentation + MuscleActiveForceRepresentation ;

% Tendon
for i = 1 : nMuscles
    for ii= i : definition
        TendonForceRepresentation(i,ii) = full(casadiFun.representationTendonForce(1,tendonLengtheningRepresentation(ii), f0m(i))) ;
    end
end

% Real
MuscleActiveForce = full(casadiFun.getMuscleActiveForce(all_states_num,muscle_tendon_parameters_num)) ;
MusclePassiveForce= full(casadiFun.getMusclePassiveForce(all_states_num,muscle_tendon_parameters_num)) ;
MuscleForce =  MuscleActiveForce + MusclePassiveForce;
TendonForce = full(casadiFun.getTendonForce(all_states_num,muscle_tendon_parameters_num)) ;



%% figures
    % muscles
figure ('Name','Mucles forces representation', 'Color',[1 1 1])
for i = 1 : nMuscles
    subplot(1,nMuscles,i)
    plot(fiberLengthRepresentation, MuscleActiveForceRepresentation(i,:) ,':k')
    hold on
    plot(fiberLengthRepresentation, MusclePassiveForceRepresentation(i,:) ,':k')
    plot(fiberLengthRepresentation, MuscleForceRepresentation(i,:) ,'k')

    plot(fiberLength(i),MuscleActiveForce(i),'or')
    plot(fiberLength(i),MusclePassiveForce(i),'or')
    plot(fiberLength(i),MuscleForce(i),'or')

    legend([{'Active muscle force [CC]'},{'Passive muscle force [CEP]'},{'Total muscle force [CC + CEP]'}])
    title(NameMuscles{i})
    xlabel('Normalized fiber lenght', 'FontWeight','bold')
    ylabel('Force (N)', 'FontWeight','bold')
    xlim([0 1.6])
    grid on
end

    % tendon
figure ('Name','Tendon forces representation', 'Color',[1 1 1])
for i = 1 : nMuscles
    subplot(1,nMuscles,i)
    plot(tendonLengtheningRepresentation + 1, TendonForceRepresentation(i,:) ,'k')
    hold on
    plot(tendonLength(i),TendonForce(i),'or')

    title(NameMuscles{i})
    xlabel('Normalized tendon lenght', 'FontWeight','bold')
    ylabel('Force (N)', 'FontWeight','bold')
    xlim([0.95 1.05])
    grid on
end

%% print in command 
for i = 1 : nMuscles
    fprintf([' states the UMT  ', NameMuscles{i}, ' : \n'])
    fprintf(['  Muscle  ', ': \n'])
    fprintf(['  -   Fiber length : ', num2str(fiberLength(i)), ' \n'])
    fprintf(['  -   CC force : ', num2str(MuscleActiveForce(i)),' Newtons', ' \n'])
    fprintf(['  -   CEP force : ', num2str(MusclePassiveForce(i)),' Newtons', ' \n'])
    fprintf(['  -   CC + CEP  force : ', num2str(MuscleForce(i)),' Newtons', ' \n'])
    fprintf(['      -->   CC force : ', num2str((MuscleActiveForce(i)/f0m(i))*100),' percent of the maximal isometric force', ' \n'])

    fprintf(' \n')
    fprintf(['  Tendon  ', ': \n'])
    fprintf(['  -   Tendon length : ', num2str(tendonLength(i)), ' \n'])
    fprintf(['  -   CES force : ', num2str(TendonForce(i)),' Newtons', ' \n'])
    fprintf(' ----------------------------------------------- \n')
end 

end
