%% forces illustration function 
% it shows force length relation ship of all muscle and tendon 
plotDefinition = 100 ; 
Length_Muscle_Relationship = linspace(0.2, 1.4, plotDefinition) ;  
Length_Tendon_Relationship = linspace(0.95, 1.05, plotDefinition) ; 

% initialisation (for perf)
Noramized_Fiber_Length_num = nan(nMuscles) ;
Noramized_Tendon_Length_num = nan(nMuscles) ;
% extract muscle parameters to do plots 
for  i =  1 : nMuscles
Noramized_Fiber_Length_num(i) = full(fiberLength_num(i) / l0m_num(i)) ; 
Noramized_Tendon_Length_num(i) = (full(tendonLengthening_num(i) )+ lst_num(i))/lst_num(i) ; 
end

% initialisation (for perf)
FMAR = nan(nMuscles) ; FMPR = nan(nMuscles) ; FMR = nan(nMuscles) ; FTR = nan(nMuscles) ;
Active_Muscle_Force_Relationship = nan(nMuscles,plotDefinition) ; 
Passive_Muscle_Force_Relationship = nan(nMuscles,plotDefinition) ; 
Muscle_Force_Relationship = nan(nMuscles,plotDefinition) ; 
Tendon_Force_Relationship = nan(nMuscles,plotDefinition) ; 

% component length to Forces 
for i =  1 : nMuscles
    Active_Muscle_Force_Relationship(i,:) = MuscleActiveForceLength_f(Length_Muscle_Relationship,f0m_num(i), a_num(i)) ; 
    Passive_Muscle_Force_Relationship(i,:) = musclePassiveForce_f(Length_Muscle_Relationship,f0m_num(i)) ; 
    Muscle_Force_Relationship(i,:) = Active_Muscle_Force_Relationship(i,:) + Passive_Muscle_Force_Relationship(i,:) ;
    Tendon_Force_Relationship(i,:) = tendonForce_f(Length_Tendon_Relationship,f0m_num(i)) ; 

    FMAR(i) = MuscleActiveForceLength_f(Noramized_Fiber_Length_num(i),f0m_num(i), a_num(i)) ; 
    FMPR(i)= musclePassiveForce_f(Noramized_Fiber_Length_num(i),f0m_num(i)); 
    FMR(i) =  FMAR(i) + FMPR(i); 
    FTR(i) = tendonForce_f(Noramized_Tendon_Length_num(i),f0m_num(i)); 

end

% plots of muscle forces 
figure ('Name','Mucles forces representation')
for i = 1 : nMuscles
subplot(1,nMuscles,i)
plot(Length_Muscle_Relationship, Active_Muscle_Force_Relationship(i,:),':k')
hold on 
plot(Length_Muscle_Relationship, Passive_Muscle_Force_Relationship(i,:),':k')
plot(Length_Muscle_Relationship, Muscle_Force_Relationship(i,:),'k')

plot(Noramized_Fiber_Length_num(i),FMAR(i),'or')
plot(Noramized_Fiber_Length_num(i),FMPR(i),'or')
plot(Noramized_Fiber_Length_num(i),FMR(i),'or')
legend([{'Active muscle force [CC]'},{'Passive muscle force [CEP]'},{'Total muscle force [CC + CEP]'}])
title([NameMuscles{i}, '- EMG : ', num2str(a_num(i))])
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
plot(Noramized_Tendon_Length_num(i), FTR(i),'or')
title([NameMuscles{i}, '- EMG : ', num2str(a_num(i))])
xlabel('Normalized lenght', 'FontWeight','bold')
ylabel('Force (N)', 'FontWeight','bold')

grid on 
end


PC_Max_Force = nan(nMuscles) ;
% muscle force in percent of max muscle force
for i = 1 : nMuscles
PC_Max_Force(i) = (FMAR(i) / f0m_num(i)) *100 ; 
end

% print in command window
fprintf('Muscle are resprectively at \n ')
for i = 1 : nMuscles
fprintf([ '   - ', NameMuscles{i}, ' : ' , num2str(PC_Max_Force(i)) , ' percent of his maximal force. \n'])
end