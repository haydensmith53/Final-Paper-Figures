%% Speed Difference Analysis - In response to reviewer comments on Smith et al. 2021 (Use through MultiSelectHayden.m)

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

for a=1:length(morphometrics.MassFromTL);
    if strcmp(INFO.whaleName,morphometrics.DeployID{a}) == 1;
        BodyMass = morphometrics.MassFromTL(a);
    end
end
clear a;
    
for a=1:length(AllData.BeatStart);
    if strcmp(INFO.whaleName,AllData.DeployID{a}) == 1;
        AllData.BodyMass{a} = BodyMass;
        AllData.BeatLength{a} = (AllData.BeatEnd(a)-AllData.BeatStart(a))/10;
        AllData.BeatStartSpd{a} = jigspeed(AllData.BeatStart(a));
        AllData.BeatEndSpd{a} = jigspeed(AllData.BeatEnd(a));
        AllData.SpdChange{a} = AllData.BeatEndSpd{a}-AllData.BeatStartSpd{a};
        AllData.SpdChngPerc{a} = abs(AllData.SpdChange{a})/AllData.BeatStartSpd{a}*100;
        AllData.DistTravel{a} = (0.5*(AllData.SpdChange{a}/AllData.BeatLength{a})*(AllData.BeatLength{a}^2))+(AllData.BeatStartSpd{a}*AllData.BeatLength{a}); % equation given by Jean
        AllData.ThrustDragOffset{a} = ((AllData.BodyMass{a}/2)*((AllData.BeatEndSpd{a}^2)-(AllData.BeatStartSpd{a}^2)))/AllData.DistTravel{a}; % equation given by Jean
        AllData.PercentOfThrust{a} = (AllData.ThrustDragOffset{a}/AllData.Thrust(a))*100;
        %plot(jigspeed(AllData.BeatStart(a):AllData.BeatEnd(a))); % plots each tailbeat as a check
        %saveas(gcf,['C:\Users\William Gough\Documents\Academic Materials\Stanford University\Github Repositories\Final-Paper-Figures\Test Beats\' INFO.whaleName ' beat#' num2str(a) '.png']);
        %clf;        
    end
end
clear a;