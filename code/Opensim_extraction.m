function [known_parameters_num,muscle_tendon_parameters_num] = Opensim_extraction()
% this function is for extract from opensim scale subjects -->
%   - known_parameters_num :
%       - segment_geometry_num
%       - muscle_origin_num
%       - muscle_insersion_num
% output type (30,1) : known_parameters_num = [segment_geometry_num; muscle_origin_num; muscle_insersion_num];
% output type (33,1) : known_parameters_num = [segment_geometry_num; muscle_origin_num; muscle_insersion_num, muscle_viaPoint];

%
%   - muscle_tendon_parameters_num :
%       - l0m_num
%       - phi0_num
%       - f0m_num
%       - lst_num
% output type (12,1) : muscle_tendon_parameters_num =  [l0m_num , phi0_num , f0m_num , lst_num ] ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% opensim file selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% name of the interst Bones
Interst_Bones = ["femur_r","tibia_r","talus_r","calcn_r","toes_r"] ;       % bones in our matlab model

% name of the interst Muscles
Interst_Muscles = ["tib_ant_r","soleus_r","lat_gas_r","med_gas_r"] ;       % muscles in our matlab model

% name of the Muscles of the model
NameModelMuscles = 'Thelen2003Muscle'  ;


CurrentFolder = cd ;
nom_fichier = [CurrentFolder, '\wholebody.osim'] ;


% ====================================
% ouverture du fichier et initialisation
% ====================================
fichier = fopen(nom_fichier,'rb');
if (fichier == -1)
    error('impossible d''ouvrir le fichier...');
end
% charger tout le fichier en entier
txt = textscan(fichier,'%s','Delimiter','\n');
fclose(fichier);
txt = char(txt{:});
N   = size(txt,1); % file length


%% define bones
nameBones = {} ; indexStartBones = [] ; indexEndBones = [] ;
nbBones = 0 ;
for i_ =  1  :  N
    if contains(string(txt(i_,:)),'<Body name')                    % Increment bone number
        nbBones = nbBones + 1 ;                                    % bone number
        indexStartBones(nbBones) = i_ ;                                     % Store the index of the description
        temp = txt(i_,:) ;                                         % Extract the char vector
        temp = regexprep(temp, '<Body name="', '');                % Remove '<Body name=' and '"' using regexprep
        temp = strrep(temp, '">', '');                             % Remove the trailing '>'
        temp = strtrim(temp);                                      % Remove leading and trailing spaces
        nameBones{nbBones} = string(temp);                         % Store the modified string
    end
    if contains(string(txt(i_,:)),'</Body>') % end of the segment desciption
        indexEndBones(nbBones) = i_ ;
    end
end
fprintf (['Bones number : ', num2str(nbBones), ' \n '])




%% define joint <CustomJoint name="
% Location of the joint in the parent body specified in the parent reference frame. Default is (0,0,0).
% <location_in_parent>
% </location_in_parent>

% extract the inforamtions of the bones
Index = find(ismember(string(nameBones),Interst_Bones));
nbJoints = 0 ; indexStartJoints = [] ; indexEndJoints = [] ; 
nameJoints = {} ; 

for i =  Index
    for ii = indexStartBones(i) : indexEndBones(i)
        if contains(string(txt(ii,:)),'<CustomJoint name="')
            nbJoints = nbJoints + 1 ;
            indexStartJoints(nbJoints) = ii ;                              % Store the index of the description
            temp = txt(ii,:) ;                                             % Extract the char vector
            temp = regexprep(temp, '<CustomJoint name="', '');             % Remove '<Body name=' and '"' using regexprep
            temp = strrep(temp, '">', '');                                 % Remove the trailing '>'
            temp = strtrim(temp);                                          % Remove leading and trailing spaces
            nameJoints{nbJoints} = string(temp);                           % Store the modified string
        end

        if contains(string(txt(ii,:)),'</Joint> ')
            indexEndJoints(nbJoints) = ii ;                                % Store the index of the end of the description
        end
    end
end

for i =  1 : nbJoints
     currentJoint = char(nameJoints{i}) ;
    for ii = indexStartJoints(i) : indexEndJoints(i)
%        string(txt(ii,:))
       
       % find location in parent 
        if contains(string(txt(ii,:)),'<location_in_parent>')
               temp = txt(ii,:) ;

                % Use regular expression to extract the three numerical
                % values of the location X Y Z
                match = regexp(temp, '<location_in_parent>([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)</location_in_parent>', 'tokens', 'once');

                if ~isempty(match)
                    joint.(currentJoint).X = str2double(match{1}) ;
                    joint.(currentJoint).Y = str2double(match{2}) ;
                    joint.(currentJoint).Z = str2double(match{3}) ;
                else
                    disp('Values not found in the string.');
                end
        end

        if contains(string(txt(ii,:)),'<parent_body>')
               temp = txt(ii,:) ;
               temp = regexprep(temp, '<parent_body>', '');                % Remove '<parent_body>' using regexprep
               temp = strrep(temp, '</parent_body>', '');                  % Remove '</parent_body>' using regexprep
               temp = strtrim(temp);                                       % Remove leading and trailing spaces
               joint.(currentJoint).parent = string(temp);                 % Store the modified string
        end
    end
end

%% define muscle </Thelen2003Muscle>
% Name and index of the desciption of the model
nameMuscles = {} ; indexStartMuscles = [] ; indexEndMuscles = [] ;
nbMuscles = 0 ;
for i_ =  1  :  N
    if contains(string(txt(i_,:)),['<',NameModelMuscles])          % Increment muscle number
        nbMuscles = nbMuscles + 1 ;                                % muscle number
        indexStartMuscles(nbMuscles) = i_ ;                        % Store the index of the description
        temp = txt(i_,:) ;                                         % Extract the char vector
        temp = regexprep(temp,['<',NameModelMuscles,' name="'],'');% Remove '<Body name=' and '"' using regexprep
        temp = strrep(temp, '">', '');                             % Remove the trailing '>'
        temp = strtrim(temp);                                      % Remove leading and trailing spaces
        nameMuscles{nbMuscles} = string(temp);                         % Store the modified string
    end
    if contains(string(txt(i_,:)),['</',NameModelMuscles]) % end of the segment desciption
        indexEndMuscles(nbMuscles) = i_ ;

    end
end

% extract the inforamtions of the muscle
Index = find(ismember(string(nameMuscles),Interst_Muscles));

% Muscle parameters
%%%%%%%%%%%%%%%%%%%%%%%%
for i =  Index
    for ii = indexStartMuscles(i) : indexEndMuscles(i)
        % tendon_slack_length
        if contains(string(txt(ii,:)),'tendon_slack_length')
            match = regexp(txt(ii,:), '<tendon_slack_length>([\d.]+)</tendon_slack_length>', 'tokens', 'once');
            if ~isempty(match)
                parameters.(char(nameMuscles{i})).tendon_slack_length = str2double(match) ; 
            else
                disp('Number not found in the string. --> tendon_slack_length \n');
            end
        end

        % optimal_fiber_length
        if contains(string(txt(ii,:)),'optimal_fiber_length')
            match = regexp(txt(ii,:), '<optimal_fiber_length>([\d.]+)</optimal_fiber_length>', 'tokens', 'once');
            if ~isempty(match)
                parameters.(char(nameMuscles{i})).optimal_fiber_length = str2double(match) ;

            else
                disp('Number not found in the string. --> optimal_fiber_length \n');
            end
        end
        %max_isometric_force
        if contains(string(txt(ii,:)),'max_isometric_force')
            match = regexp(txt(ii,:), '<max_isometric_force>([\d.]+)</max_isometric_force>', 'tokens', 'once');
            if ~isempty(match)
                parameters.(char(nameMuscles{i})).max_isometric_force = str2double(match) ;
            else
                disp('Number not found in the string. --> max_isometric_force \n');
            end
        end

        % pennation_angle_at_optimal
        if contains(string(txt(ii,:)),'pennation_angle_at_optimal')
            match = regexp(txt(ii,:), '<pennation_angle_at_optimal>([\d.]+)</pennation_angle_at_optimal>', 'tokens', 'once');
            if ~isempty(match)
                parameters.(char(nameMuscles{i})).pennation_angle_at_optimal = str2double(match) ;
            else
                disp('Number not found in the string. --> pennation_angle_at_optimal \n');
            end
        end
    end
end

% Muscle insertions
%%%%%%%%%%%%%%%%%%%%%%%%
nbPathPoint = 0 ;

% PathPoint location in the opensim file
for i =  Index
    for ii = indexStartMuscles(i) : indexEndMuscles(i)
        if contains(string(txt(ii,:)),'<PathPoint name=')
            % PathPoint index
            nbPathPoint = nbPathPoint + 1 ;                                % muscle number
            indexStartPathPoint(nbPathPoint) = ii;                         % Store the index of the description
        end

        if contains(string(txt(ii,:)),'</PathPoint>')
            % index
            indexEndPathPoint(nbPathPoint) = ii;                           % Store the index of the description
        end
    end
    
    % PathPoint information
    insertion = [];
    for i = 1: nbPathPoint
        temp = txt(indexStartPathPoint(i),:) ;
        temp = regexprep(temp,'<PathPoint name="','');                     % Remove '<Body name=' and '"' using regexprep
        temp = strrep(temp, '">', '')    ;                                 % Remove the trailing '>'
        insertionName = string(strtrim(temp));                             % Remove leading and trailing spaces
        insertionName = strrep(insertionName, '-', '_') ;

        for ii = indexStartPathPoint(i) : indexEndPathPoint(i)
            % PathPoint location <location>
            if contains(string(txt(ii,:)),'<location>')
                temp = txt(ii,:) ;

                % Use regular expression to extract the three numerical
                % values of the location X Y Z
                match = regexp(temp, '<location>\s+([\d.-]+)\s+([\d.-]+)\s+([\d.-]+)</location>', 'tokens', 'once') ;

                if ~isempty(match)
                    insertion.(insertionName).X= str2double(match{1}) ;
                    insertion.(insertionName).Y= str2double(match{2}) ;
                    insertion.(insertionName).Z= str2double(match{3}) ;
                else
                    disp('Values not found in the string.');
                end
            end

            % PathPoint body <body>
            if contains(string(txt(ii,:)),'<body>')
                temp = txt(ii,:) ; 
                temp = regexprep(temp,'<body>','');                      % Remove '<Body name=' and '"' using regexprep
                temp = strrep(temp, '</body>', '')    ;                    % Remove the trailing '>'
                insertion.(insertionName).parent=  string(strtrim(temp));  % Remove leading and trailing spaces
            end
        end
    end
end

    %% output of our model
    try
            % gast : sum of lat_gas_r and med_gas_r
                %insertion 
            insertion.gast_P1.X = mean([insertion.med_gas_r_P1.X, insertion.lat_gas_r_P1.X]) ;
            insertion.gast_P1.Y = mean([insertion.med_gas_r_P1.Y, insertion.lat_gas_r_P1.Y]) ;
            insertion.gast_P1.Z = mean([insertion.med_gas_r_P1.Z, insertion.lat_gas_r_P1.Z]) ; 
            insertion.gast_P1.parent = insertion.med_gas_r_P1.parent ; 

            insertion.gast_P3.X = mean([insertion.med_gas_r_P3.X, insertion.lat_gas_r_P3.X]) ;
            insertion.gast_P3.Y = mean([insertion.med_gas_r_P3.Y, insertion.lat_gas_r_P3.Y]) ;
            insertion.gast_P3.Z = mean([insertion.med_gas_r_P3.Z, insertion.lat_gas_r_P3.Z]) ; 
            insertion.gast_P3.parent = insertion.med_gas_r_P3.parent ; 

                %parameters : F0m --> sum ; L0m --> mean ; phi0 --> mean , tsl --> mean 
            
           parameters.gast.optimal_fiber_length = mean([parameters.med_gas_r.optimal_fiber_length, ... 
               parameters.lat_gas_r.optimal_fiber_length]) ; 
           parameters.gast.pennation_angle_at_optimal = mean([parameters.med_gas_r.pennation_angle_at_optimal, ...
               parameters.lat_gas_r.pennation_angle_at_optimal]) ;
           parameters.gast.max_isometric_force = sum([parameters.med_gas_r.max_isometric_force, ... 
               parameters.lat_gas_r.max_isometric_force]) ; 
           parameters.gast.tendon_slack_length = mean([parameters.med_gas_r.tendon_slack_length, ... 
               parameters.lat_gas_r.tendon_slack_length]) ; 

        % known_parameters_num = [segment_geometry_num; muscle_origin; muscle_insersion];
            %segment_geometry_num = [thigh_0, thigh_1, thigh_2, leg_0, leg_1, leg_2, talus_0, talus_1, talus_2, foot_0, foot_1, foot_2]
            segment_geometry_num = [joint.knee_r.X, joint.knee_r.Y, 0,...
                joint.ankle_r.X, joint.ankle_r.Y, 0,...
                joint.subtalar_r.X, joint.subtalar_r.Y, 0,...
                joint.mtp_r.X, joint.mtp_r.Y, 0] ; % OK 

%             muscle_origin = [Local_Origin_tibialis_anterior_0, Local_Origin_tibialis_anterior_1, Local_Origin_tibialis_anterior_2, Local_Origin_soleus_0, Local_Origin_soleus_1, Local_Origin_soleus_2, Local_Origin_gastrocnemius_0, Local_Origin_gastrocnemius_1, Local_Origin_gastrocnemius_2, Local_Insertion_tibialis_anterior_0, Local_Insertion_tibialis_anterior_1, Local_Insertion_tibialis_anterior_2, Local_Insertion_soleus_0, Local_Insertion_soleus_1, Local_Insertion_soleus_2, Local_Insertion_gastrocnemius_0, Local_Insertion_gastrocnemius_1, Local_Insertion_gastrocnemius_2] ;
            muscle_origin = [insertion.tib_ant_r_P1.X, insertion.tib_ant_r_P1.Y, 0, ...
                insertion.soleus_r_P1.X, insertion.soleus_r_P1.Y, 0, ...
                insertion.gast_P1.X, insertion.gast_P1.Y, 0] ; % OK

            muscle_insersion = [insertion.tib_ant_r_P3.X, insertion.tib_ant_r_P3.Y, 0, ...
                insertion.soleus_r_P2.X, insertion.soleus_r_P2.Y, 0, ...
                insertion.gast_P3.X, insertion.gast_P3.Y, 0] ; 

            muscle_viaPoint = [insertion.tib_ant_r_P2.X, insertion.tib_ant_r_P2.Y, 0] ;

known_parameters_num = [segment_geometry_num, muscle_origin, muscle_insersion, muscle_viaPoint ] ; 

        % muscle_tendon_parameters_num =  [l0m_num , phi0_num , f0m_num , lst_num ] ;
            l0m_num = [parameters.tib_ant_r.optimal_fiber_length,...
                parameters.soleus_r.optimal_fiber_length,...
                parameters.gast.optimal_fiber_length ] ;  
            phi0_num = [parameters.tib_ant_r.pennation_angle_at_optimal,...
                parameters.soleus_r.pennation_angle_at_optimal,...
                parameters.gast.pennation_angle_at_optimal ] ;  
            f0m_num = [parameters.tib_ant_r.max_isometric_force,...
                parameters.soleus_r.max_isometric_force,...
                parameters.gast.max_isometric_force ] ; 
            lst_num =  [parameters.tib_ant_r.tendon_slack_length,...
                parameters.soleus_r.tendon_slack_length,...
                parameters.gast.tendon_slack_length ] ; 

muscle_tendon_parameters_num =  [l0m_num , phi0_num , f0m_num , lst_num ] ;  

fprintf(['Optimal fiber length (l0m) [Tibialis Anterior, Soleus, Gastrocnemius]','\n '])
fprintf([num2str(l0m_num),'\n'])

fprintf(['Pennation angle at optimal fiber length (phi0) [Tibialis Anterior, Soleus, Gastrocnemius] ','\n '])
fprintf([num2str(phi0_num),'\n'])

fprintf(['Maximal fiber force (f0m) [Tibialis Anterior, Soleus, Gastrocnemius]','\n '])
fprintf([num2str(f0m_num),'\n'])

fprintf(['tendon slack length (lst) [Tibialis Anterior, Soleus, Gastrocnemius]','\n '])
fprintf([num2str(lst_num),'\n'])

    catch
        known_parameters_num = zeros(1,33) ; 
        muscle_tendon_parameters_num = zeros(1,33) ; 
        fprintf('error in parameter extraction : function Opensim_extraction /n')
    end


end