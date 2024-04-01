%% protocol

%% Passives conditions (mobilisations) 
PassivesConditions = Data(all(Data(:,4:6) == 0,2),:) ; % Neuromuscular activation = 0
PassivesConditions = [PassivesConditions(PassivesConditions(:,2) == (0/180)*pi, :) ; ... % knee angle = 0 deg
    PassivesConditions(PassivesConditions(:,2) == (60/180)*pi, :)] ; % knee angle = 60 deg

%% Actives conditions
%   iMVC (a =1 )
%   q1 : 0 and 80 
%   q2 : [-20, - 10, 0 , 10 , 20 ]

ActivesConditions = Data(any(Data(:,4:6) == 1,2),:) ; % Neuromuscular activation = 0
TibialisAnteriorActivesConditions = ActivesConditions(ActivesConditions(:,4) == 1,:) ; % Neuromuscular activation = 1
SoleusGacstrocnemiusActivesConditions = ActivesConditions(any(ActivesConditions(:,5:6) == 1,2),:) ; % Neuromuscular activation = 1

TibialisAnteriorActivesConditions = [TibialisAnteriorActivesConditions(TibialisAnteriorActivesConditions(:,2) == (0/180)*pi, :)] ; % knee angle = 60 deg

SoleusGacstrocnemiusActivesConditions = [SoleusGacstrocnemiusActivesConditions(SoleusGacstrocnemiusActivesConditions(:,2) == (0/180)*pi, :) ; ... % knee angle = 0 deg
    SoleusGacstrocnemiusActivesConditions(SoleusGacstrocnemiusActivesConditions(:,2) == (60/180)*pi, :)] ; % knee angle = 60 deg

qankle = ([-20, - 10, 0 , 10 , 20 ]/ 180)*pi ; 
SelectedActivesConditions = [] ; 

for i =  1 : size(TibialisAnteriorActivesConditions,1)
    if any(TibialisAnteriorActivesConditions(i,3) == qankle)
        SelectedActivesConditions = [SelectedActivesConditions ; TibialisAnteriorActivesConditions(i,:)];
    end 
end

for i =  1 : size(SoleusGacstrocnemiusActivesConditions,1)
    if any(SoleusGacstrocnemiusActivesConditions(i,3) == qankle)
        SelectedActivesConditions = [SelectedActivesConditions ; SoleusGacstrocnemiusActivesConditions(i,:)];
    end 
end
nActivesTrials = size(SelectedActivesConditions,1); 

SelectedConditions = [PassivesConditions; ... 
    SelectedActivesConditions] ; 



%% sortir les tendons forces pour chaque Q1 et Q2 