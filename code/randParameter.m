function [pRand] = randParameter(pOrigin,Noise)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1
    % Noise
    %%%%%%%%%%%%%%%%%%%%%
    l0mNoise = .5;
    f0mNoise = .5;
    tslNoise = .5;

else
    l0mNoise = Noise;
    f0mNoise = Noise;
    tslNoise = Noise;
end


nMuscle = 3;
% parameters
%%%%%%%%%%%%%%%%%%%%%
l0m = pOrigin(1:nMuscle);
% alpha = pOrigin(nMuscle+1:nMuscle*2);
f0m = pOrigin((nMuscle*2)+1:nMuscle*3);
tsl = pOrigin((nMuscle*3)+1:nMuscle*4);


% range in random values
%%%%%%%%%%%%%%%%%%%%%%%%
l0mRange = [l0m - l0m*l0mNoise; l0m + l0m*l0mNoise];                       % in meter
angleRange = ones(2,nMuscle) .* (([10; 60])./180).*pi;                      % in radian
f0mRange = [f0m - f0m*f0mNoise; f0m + f0m*f0mNoise];                       % in Newton
tslRange = [tsl - tsl*tslNoise; tsl + tsl*tslNoise];                       % in meter

Range = [l0mRange, angleRange, f0mRange,tslRange];
randKey = rand(1,size(pOrigin,2));

pRand = NaN(1,size(pOrigin,2));
for i= 1 : size(Range,2)
    pRand(i) = randInRange(Range(:,i),randKey(i));
end

end 

function val = randInRange(range,randValue)

    val = range(1) + (range(2) - range(1)) * randValue;
end