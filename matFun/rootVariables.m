function [rooted_variables,temp_residuals] = rootVariables(a,lmtu,p,casadiFun,toleratedError,nTry,unknown)
if nargin == 6
    unknown  = bestStart(a,lmtu,p);
end

known_num = [a,lmtu,p];

currentTry = 0;
rooted_variables = NaN(1,5);
err_mtu_length = NaN(1);
err_forces = NaN(1);

temp_err = 100;
temp_sol = NaN(1);
temp_err_mtu_length = NaN(1);
temp_pc_err_forces = NaN(1);
temp_residuals = NaN(1);
success = false;

while currentTry < nTry
    % to solve 
    residuals = full(casadiFun.equilibriumErrorSingleMuscle(unknown',known_num)); % residuals
    rooted_variables= full(casadiFun.equilibrateMuscleTendonSingleMuscle(unknown,known_num));

    rooted_mtu_length = cos (rooted_variables(5)) * rooted_variables(4) +  rooted_variables(3) + p(4);
    err_mtu_length = abs(rooted_mtu_length - known_num(2));
    err_forces = abs(rooted_variables(1) - rooted_variables(2));

    pc_err_mtu_length = (err_mtu_length/known_num(2))*100;
    pc_err_forces = (err_forces/mean(rooted_variables(1:2)))*100;

    if temp_err > mean(abs([pc_err_mtu_length, pc_err_forces]))
        temp_err = mean(abs([pc_err_mtu_length, pc_err_forces]));
        temp_err_mtu_length = pc_err_mtu_length;
        temp_pc_err_forces = pc_err_forces;
        temp_residuals = residuals;
        temp_sol = rooted_variables;
    end


    if err_mtu_length<toleratedError && err_forces<toleratedError          %all(abs(residuals(3:end))<toleratedError) % find a solution 
        rooted_variables= full(casadiFun.equilibrateMuscleTendonSingleMuscle(unknown,known_num));
        disp(['number of try ', num2str(nTry)])
        disp(['sol find: residuals ', num2str(toleratedError)])
        disp(['residuals ', num2str(residuals')])

        fprintf('\n\n ')
        disp('  Error in equilibrium : ')
        disp(['error in length mtu : ', num2str(err_mtu_length)])
        disp(['error in forces : ', num2str(err_forces)])
        success = true;

        break
    end 

    fprintf('\n\n ')
    disp('  Error in equilibrium : ')
    disp(['error in length mtu (in m) : ', num2str(err_mtu_length)])
    disp(['error in length mtu (in %): ', num2str(pc_err_mtu_length)])

    disp(['error in forces (in N) : ', num2str(err_forces)])    
    disp(['error in forces (in %): ', num2str(pc_err_forces)])

    % pause

    % clear 
    rooted_variables = NaN(1,5);
    currentTry = currentTry + 1;
    unknown  = bestStart(a,lmtu,p);

end

%% output 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if success == true
    TendonForce = rooted_variables(1);
    MuscleForce = rooted_variables(2);
    tendon_lengthening_rooted = rooted_variables(3);
    fiber_length_rooted = rooted_variables(4);
    pennation_angle_rooted = rooted_variables(5);
else
    TendonForce = temp_sol(1);
    MuscleForce = temp_sol(2);
    tendon_lengthening_rooted = temp_sol(3);
    fiber_length_rooted = temp_sol(4);
    pennation_angle_rooted = temp_sol(5);
end
pennation_deg = (pennation_angle_rooted/pi)*180;
tendon_length_rooted = p(4) + tendon_lengthening_rooted;

fprintf('\n\n ')
fprintf('Results \n')
fprintf('===================================================================')
disp('  about algo : ')
disp(['Tolerated error  : ', num2str(toleratedError)])
disp(['Number of test  : ', num2str(currentTry), ' out of ', num2str(nTry)])
disp(['Mean error : ', num2str(temp_err)])
disp(['state : ', num2str(success)])



fprintf('\n\n ')
disp('  Error in equilibrium : ')
disp(['error in length mtu (in m) : ', num2str(err_mtu_length)])
disp(['error in length mtu (in %): ', num2str(temp_err_mtu_length)])

disp(['error in forces (in N) : ', num2str(err_forces)])
disp(['error in forces (in %): ', num2str(temp_pc_err_forces)])

fprintf('\n\n')
disp(['rooted fiber length: ', num2str(fiber_length_rooted),' m'])
disp(['rooted tendon length: ', num2str(tendon_length_rooted),' m'])
disp(['rooted pennation angle: ', num2str(pennation_deg),' deg'])
disp(['rooted tendon forces: ', num2str(TendonForce),' N'])
disp(['rooted muscle forces: ', num2str(MuscleForce),' N'])




end 