function [pRand] = randParameter(pOrigin,Noise)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1
    % Noise
    %%%%%%%%%%%%%%%%%%%%%
    l0mNoise = .5;
    f0mNoise = .5;
    phi0Noise = .5;
    tslNoise = .5;

else
    if length(Noise) == 1
        l0mNoise = Noise;
        f0mNoise = Noise;
        tslNoise = Noise;
        phi0Noise = Noise;
    elseif length(Noise) == 4
        l0mNoise = Noise(1);
        f0mNoise = Noise(2);
        tslNoise = Noise(3);
        phi0Noise = Noise(4);
    else
        error('length of Noise must be (1x1) or (1x4)')
    end
end


nMuscle = 3;
% parameters
%%%%%%%%%%%%%%%%%%%%%
l0m = pOrigin(1:nMuscle);
phi0 = pOrigin(nMuscle+1:nMuscle*2);
f0m = pOrigin((nMuscle*2)+1:nMuscle*3);
tsl = pOrigin((nMuscle*3)+1:nMuscle*4);


% range in random values
%%%%%%%%%%%%%%%%%%%%%%%%
l0mRange = [l0m - l0m*l0mNoise; l0m + l0m*l0mNoise];                       % in meter
phi0Range = [phi0 - phi0*phi0Noise; phi0 + phi0*phi0Noise];                     % in radian
f0mRange = [f0m - f0m*f0mNoise; f0m + f0m*f0mNoise];                       % in Newton
tslRange = [tsl - tsl*tslNoise; tsl + tsl*tslNoise];                       % in meter

% add border to the noise to avoid impossible value (e.g. negative forces,
% negative angle, negative length)
%%%%%%%%%%%%%%%%%%%%%%%%

l0mRange(l0mRange<=0) = 0.0001;                                            % in meter
phi0Range(phi0Range<=(deg2rad(10))) = deg2rad(10);                         % in radian
phi0Range(phi0Range>=(deg2rad(80))) = deg2rad(80);                         % in radian
f0mRange(f0mRange<=0) = 1;                                                 % in Newton
tslRange(tslRange<=0) = 0.0001; 

% create random values in ranges
%%%%%%%%%%%%%%%%%%%%%%%%%
Range = [l0mRange, phi0Range, f0mRange,tslRange];
randKey = rand(1,size(pOrigin,2));

pRand = NaN(1,size(pOrigin,2));
for i= 1 : size(Range,2)
    pRand(i) = randInRange(Range(:,i),randKey(i));
end
end 