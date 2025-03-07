function plotTendonForcelength(vizualizationFun,tendon_slack_length,maximalIsometricForce,tendon_current_length)
% check the input
%%%%%%%%%%%%%%%%%%
if nargin == 1
    tendon_slack_length = 1;
    maximalIsometricForce = 1;
    tendon_current_length = NaN(1);
elseif nargin == 3
    tendon_current_length = NaN(1);
end

% about lenght
tendon_normalized_length = linspace(0.95,1.05,1000);
tendon_length = tendon_normalized_length * tendon_slack_length;
tendon_current_normalized_length = tendon_current_length/tendon_slack_length;

% about force
tendon_normalized_force = full(vizualizationFun.getNormalizedTendonForce(tendon_length,tendon_slack_length));
tendon_force = tendon_normalized_force * maximalIsometricForce;
tendon_current_normalized_force = full(vizualizationFun.getNormalizedTendonForce(tendon_current_length,tendon_slack_length));
tendon_current_force = tendon_current_normalized_force * maximalIsometricForce;

%% plot 
figure("Name","Tendon force-length relationship")
subplot(1,2,1)
hold on 
title('Normalized force-lenght relationship')
plot(tendon_normalized_length,tendon_normalized_force,'k')
plot(tendon_current_normalized_length,tendon_current_normalized_force,'or')
xlabel('Normalized tendon length','FontWeight','bold')
ylabel('Normalized tendon force','FontWeight','bold')
xlim([min(tendon_normalized_length) max(tendon_normalized_length)])
ylim([min(tendon_normalized_force) max(tendon_normalized_force)])
legend('force-length','current value')
hold off

subplot(1,2,2)
hold on 
title('force-lenght relationship')
plot(tendon_length,tendon_force,'k')
plot(tendon_current_length,tendon_current_force,'or')
xlabel('Tendon length (m)','FontWeight','bold')
ylabel('Tendon force (N)','FontWeight','bold')
xlim([min(tendon_length) max(tendon_length)])
ylim([min(tendon_force) max(tendon_force)])
hold off

end 