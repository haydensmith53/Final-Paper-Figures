%% Speed Difference Analysis - In response to reviewer comments on Smith et al. 2021

beforewhale = find(tagon,1); % finds the index of tag placement
afterwhale = find(tagon);
afterwhale = afterwhale(end); % finds the index where the tag falls off the whale

jigspeed = speed.JJ(beforewhale:end,:);
NaNs = find(isnan(jigspeed)==1); 
jigspeed(NaNs,:) = []; % throws out NaNs

fs = 10; % frame rate
N = 35; % order
Fpasslow = 0.09; % lowpass filter passband frequency
Fstoplow = 0.2; % lowpass filter stopband frequency
Wpasslow = 0.4; % passband weight
Wstoplow = 0.4; % stopband weight
denslow  = 60; % density factor
blow  = firpm(N, [0 Fpasslow Fstoplow fs/2]/(fs/2), [1 1 0 0], [Wpasslow Wstoplow], ...
       {denslow});
Hdlow = dfilt.dffir(blow);
jigspeed(:,2) = filtfilt(blow,1,jigspeed); % lowpass filters the jiggle speed

for a=1:length(AllData.BeatStart);
    if strcmp(INFO.whaleName,AllData.DeployID{a}) == 1;
        AllData.SpdChange{a} = jigspeed(AllData.BeatEnd(a)) - jigspeed(AllData.BeatStart(a));
    end
end
clear a;