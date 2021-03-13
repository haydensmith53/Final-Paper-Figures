%%  Large Whale Swimming Kinematics / Hydrodynamics Project
%   Hayden Smith & William Gough (2021) - updated from previous version
%   written by William Gough and Max Czapanskiy (2019).
%   The goal of this project is to determine the swimming kinematics and
%   subsequent hydrodynamics for large rorqual whales.
%clear all;
%[filename, fileloc] = uigetfile('MultiSelect','off'); % import the PRH file here
%load([fileloc filename]);

morphometrics = readtable('C:\Users\William Gough\Documents\Academic Materials\Stanford University\Github Repositories\Final-Paper-Figures\Finalized Data Sheet For Hayden.xlsx'); % read in all of the morphometric data

for a=1:length(morphometrics.DeployID); % identify the whale from the list of possible individuals/deployments in the morphometric data list
    if strcmp(INFO.whaleName,morphometrics.DeployID(a));
        TL = morphometrics.TotLength(a); % total length in meters
        mass = morphometrics.MassFromTL(a); % mass based on Shirel's allometric curves
        Sa = morphometrics.SurfArea(a); % wetted surface area in meters squared
        maxdiam = morphometrics.MaxDiam(a); % widest diameter of the animal in meters
        Fa = morphometrics.FlukeArea(a); % planform area of flukes in meters squared
        C = morphometrics.ChordLength(a); % chord length of fluke in meters
        Species = morphometrics.Species(a); % species of the current deployment
        break
    else
        continue
    end
end
clear a;

beforewhale = find(tagon,1); % finds the index of tag placement
afterwhale = find(tagon);
afterwhale = afterwhale(end); % finds the index where the tag falls off the whale

offsetzero = 0;
ptrunc = p((beforewhale):end,1); % the truncated depth using tag on/off indices
pitchtrunc = pitch((beforewhale):end,1); % the truncated pitch

if strcmp(INFO.whaleName,'bp170907-41b') || strcmp(INFO.whaleName,'bw180904-48') || strcmp(INFO.whaleName,'mn180227-41')
    yaxisgyro = Aw((beforewhale):end,1); % the truncated y-axis accelerometer signal
else
    yaxisgyro = Gw((beforewhale):end,2); % the truncated y-axis gyroscope signal
end

NaNs = find(isnan(yaxisgyro)==1); 
yaxisgyro(NaNs,:) = []; % throws out NaNs
ptrunc(NaNs,:) = [];
pitchtrunc(NaNs,:) = [];

deptharray = [];
timeabovezeroavg = [];
zeroavg = mean(yaxisgyro(100:(end-100))); % the mean gyroscope value to be used as a crossing for our "zero-crossings"
flukingstartpoints = [];

fluklengththresh = 100; % threshold for flukebeat length (i.e. the length of a flukebeat must be less than this many indices)
flukmagnitudethresh = 1; % threshold for flukebeat magnitude (i.e. the flukebeat must be at least this high on both the upstroke and downstroke)
flukclaritythresh = 15; % threshold for flukebeat clarity/symmetry (i.e. the upstroke and downstroke must be within this multiplicative value of one another)
flukheightthresh = zeroavg*2; % threshold for finding peaks in a flukebeat
    
for a=1:length(yaxisgyro); % determines periods of depth above a depth threshold (2m) 
    if ptrunc(a)<2;
        deptharray(a,1) = 0;
    else
        deptharray(a,1) = 1;
    end
end
clear a;

fs = 10; % frame rate
N = 10; % order
Fpasslow = 0.44; % lowpass filter passband frequency
Fstoplow = 1.0; % lowpass filter stopband frequency
Wpasslow = 0.1; % passband weight
Wstoplow = 0.1; % stopband weight
denslow  = 20; % density factor
blow  = firpm(N, [0 Fpasslow Fstoplow fs/2]/(fs/2), [1 1 0 0], [Wpasslow Wstoplow], ...
       {denslow});
Hdlow = dfilt.dffir(blow);
yaxisgyro(:,2) = filtfilt(blow,1,yaxisgyro); % lowpass filters the dataset

if strcmp(INFO.whaleName,'bp170907-41b') || strcmp(INFO.whaleName,'bw180904-48') || strcmp(INFO.whaleName,'mn180227-41')
    Fstophigh = 0.072;    % Highpass Filter Stopband Frequency
    Fpasshigh = 0.1368;   % Highpass Filter Passband Frequency
    Dstophigh = 0.0001;          % Stopband Attenuation
    Dpasshigh = 0.057501127785;  % Passband Ripple
    denshigh  = 20;              % Density Factor
    [N, Fo, Ao, W] = firpmord([Fstophigh, Fpasshigh]/(fs/2), [0 1], [Dstophigh, Dpasshigh]);
    bhigh  = firpm(N, Fo, Ao, W, {denshigh});
    Hdhigh = dfilt.dffir(bhigh);
    yaxisgyro(:,2) = filtfilt(bhigh,1,yaxisgyro(:,2));
end

for a=1:length(yaxisgyro(:,2)); % creates an array for segments when the gyroscope data is above the average of the dataset (proxy for fluking - "zero crossing")
    if (yaxisgyro(a,2)) >= (zeroavg); % change this between 0 and (zeroavg*x)
        timeabovezeroavg(a,1) = 1;
    else
        timeabovezeroavg(a,1) = 0;
    end
end
clear a;

flukingstartpoints(1) = 0;
lookahead = 3;

for a=2:length(yaxisgyro(:,2))-lookahead; % finds the start of individual flukebeats
    if timeabovezeroavg(a) == 1 && timeabovezeroavg(a-1) == 0 && timeabovezeroavg(a+lookahead) == 1; % tests if there is a zero crossing
        flukingstartpoints(a,1) = 1;
    else
        flukingstartpoints(a,1) = 0;
    end
end
clear a;

for a=1:lookahead % adds zeros onto the end of the flukebeat list
    flukingstartpoints(end+1) = 0;
end
clear a;

flukebeatstart = [];
flukebeatend = [];
useableflukebeats = [];

for a=1:length(flukingstartpoints); % creates a list of starting and ending indices for each flukebeat found in the previous step
    if flukingstartpoints(a) == 0;
        continue
    else
        flukebeatstart(end+1,1) = a;
        for b=(a+1):length(flukingstartpoints);
            if flukingstartpoints(b) == 0;
                continue
            else
                flukebeatend(end+1,1) = b;
                break
            end
        end
    end
end
clear a b;

flukebeatend(end+1) = length(flukingstartpoints);
useableflukebeats = [flukebeatstart flukebeatend];

for a=1:length(useableflukebeats(:,1)); % tests whether a particular flukebeat is less than a threshold length (indices). If it is, a trapezoidal integration is run for the flukebeat
    if minus((useableflukebeats(a,2)),(useableflukebeats(a,1))) >= fluklengththresh;
        continue
    else
        tempforint = yaxisgyro(((useableflukebeats(a,1)):(useableflukebeats(a,2))),2);
        useableflukebeats(a,3) = trapz(tempforint); % runs the trapezoidal integration and adds it to the third colum of the flukebeat list
        useableflukebeats(a,4) = trapz(abs(tempforint)); % runs another trapezoidal integration to determine approximate magnitude of the flukebeat
        [pks,locs] = findpeaks(yaxisgyro(useableflukebeats(a,1):useableflukebeats(a,2),2));
        [pksabs,locsabs] = findpeaks(abs(yaxisgyro(useableflukebeats(a,1):useableflukebeats(a,2),2)));
        highpks = find(pksabs > (flukheightthresh));
        if useableflukebeats(a,3) <= flukclaritythresh && ... % tests to make sure the flukebeat isn't too top-heavy
           useableflukebeats(a,3) >= -flukclaritythresh && ... % tests to make sure the flukebeat isn't too bottom-heavy
           useableflukebeats(a,4) >= flukmagnitudethresh && ... % tests to make sure the flukebeat magnitude is above a threshold
           length(locs) == 1 && ... % tests to make sure there is only one positive peak in the flukebeat
           length(locsabs) == 2 && ... % tests to make sure there are exactly two peaks (one positive, one negative) in a flukebeat
           isempty(highpks) == 0; % tests to make sure the flukebeat is at least as high as a threshold in the positive direction
            continue
        else
            useableflukebeats(a,3:4) = 0;
        end
    end
end
clear a;

indices = find(useableflukebeats(:,3)==0); % removes flukebeats that are not within desired boundaries
useableflukebeats(indices,:) = [];

isfluking = NaN(length(yaxisgyro(:,1)),1);

for flukephase = useableflukebeats.' % creates a list of indices when the whale is (1) and is not (0) fluking
    isfluking(flukephase(1):flukephase(2)) = 1;
end

yaxisgyro = [yaxisgyro isfluking]; % appends the list of fluking times onto the gyroscope data
yaxisgyro = [yaxisgyro abs(yaxisgyro(:,2))];

useableflukebeats(:,5) = minus((useableflukebeats(:,2)),(useableflukebeats(:,1))); % calculate the frequency of each flukebeat
useableflukebeats(:,6) = useableflukebeats(:,5) ./ 10;
useableflukebeats(:,7) = 1 ./ useableflukebeats(:,6); % this column is the actual oscillatory frequency of the flukebeat
dsf = mean(useableflukebeats(:,7)); % the average oscillatory frequency for the whale (all flukebeats - deep and shallow combined)

BeatStart = useableflukebeats(:,1); % flukebeat start indices (from beforewhale, not beginning of deployment)
BeatEnd = useableflukebeats(:,2); % flukebeat end indices (from beforewhale, not beginning of deployment)
OsFreq = useableflukebeats(:,7); % flukebeat oscillatory frequencies

useableflukebeatdeep = table(BeatStart,BeatEnd,OsFreq); % new table that includes flukebeat start/end indices and oscillatory frequencies

useableflukebeatdeep.BeatStartFBD = useableflukebeatdeep.BeatStart+beforewhale; % flukebeat start index from deployment beginning
useableflukebeatdeep.BeatEndFBD = useableflukebeatdeep.BeatEnd+beforewhale; % flukebeat end index from deployment beginning

for a=1:length(useableflukebeatdeep.BeatStart); % adds average depth and pitch values and a determinant of deep/shallow (deptharray)
    useableflukebeatdeep.AvgDepth(a) = mean(ptrunc((useableflukebeatdeep.BeatStart(a)):(useableflukebeatdeep.BeatEnd(a))));
    useableflukebeatdeep.AvgPitch(a) = mean(pitchtrunc((useableflukebeatdeep.BeatStart(a)):(useableflukebeatdeep.BeatEnd(a)))*180/pi);
    useableflukebeatdeep.DeepOrSurface(a) = deptharray(useableflukebeatdeep.BeatStart(a));
end
clear a;

useableflukebeatdeep(useableflukebeatdeep.DeepOrSurface == 0,:) = []; % removes shallow flukebeats
useableflukebeatdeep.DeepOrSurface = []; % removes the .DeepOrSurface column (not necessary after this)

isflukingdeep = NaN(length(yaxisgyro(:,1)),1); % recreates isflukingdeep but for deep flukebeats only

for a=1:length(useableflukebeatdeep.BeatStart); % finalizes isflukingdeep using the values from useableflukebeatdeep for flukebeat index locations
    isflukingdeep(useableflukebeatdeep.BeatStart(a):useableflukebeatdeep.BeatEnd(a)) = zeroavg;
end

yaxisgyro = [yaxisgyro isflukingdeep]; % appends the list of deep fluking times onto the gyroscope data
dsfdeep = mean(useableflukebeatdeep.OsFreq); % the average oscillatory frequency for the whale (only deep flukebeats)


disp('Dominant Stroking Frequency (Below X Meters): ');
disp(dsfdeep);
disp('Number of Tailbeats Used (Below X Meters): ');
disp(length(useableflukebeatdeep.OsFreq));


%% Find Speed Measures For Each Tailbeat
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

for a=1:length(useableflukebeatdeep.BeatStart);
    if useableflukebeatdeep.BeatEnd(a) >= length(jigspeed(:,2));
        useableflukebeatdeep((a:end),:) = [];
        break
    else
        useableflukebeatdeep.AvgSpeeds(a) = mean(jigspeed((useableflukebeatdeep.BeatStart(a)):(useableflukebeatdeep.BeatEnd(a)),2));
        useableflukebeatdeep.MedianSpeeds(a) = median(jigspeed((useableflukebeatdeep.BeatStart(a)):(useableflukebeatdeep.BeatEnd(a)),2));
    end
end
clear a;

speedmin = min(speed.JJ(beforewhale:afterwhale,1));
speedthresh = speedmin+nanstd(speed.JJ(beforewhale:afterwhale,1));

highspeedbeats = useableflukebeatdeep;
highspeedbeats(highspeedbeats.AvgSpeeds < speedthresh,:) = []; % removes shallow flukebeats

percentile25 = prctile(useableflukebeatdeep.AvgSpeeds,25);
median = prctile(useableflukebeatdeep.AvgSpeeds,50);
percentile75 = prctile(useableflukebeatdeep.AvgSpeeds,75);

if strcmp(INFO.whaleName,'bb180304-45')
    useableflukebeatdeep = sortrows(useableflukebeatdeep,'AvgSpeeds','ascend');
    useableflukebeatdeep(577:1847,:) = [];
    useableflukebeatdeep = sortrows(useableflukebeatdeep,'BeatStart','ascend');
end
  
        


%% Import Lunge File and Determine Routine Vs. Lunge-Associated Effort Tailbeats
lungedir = 'D:\REU Hayden\prhlunges\'; % this is the directory with all lunge files

try
    D = dir([lungedir '*' INFO.whaleName '*']); % this section loads the desired lunge file
    D.folder = lungedir;
    load([D.folder '\' D.name]);
    numlunges = length(LungeI);
    disp([num2str(numlunges) ' lunges in deployment ' INFO.whaleName]);
catch
end

if exist('LungeI','var') == 0;
    LungeI = NaN(1000,1);
end

thislunge = 1; % this section will determine whether each tailbeat is associated with a lunge (within 10 seconds) or routine swimming
thislungeI = LungeI(thislunge);
for a = 1:length(useableflukebeatdeep.BeatStartFBD);
    if useableflukebeatdeep.BeatStartFBD(a) >= thislungeI - 100 && useableflukebeatdeep.BeatEndFBD(a) <= thislungeI;
        useableflukebeatdeep.MaxOrNormal{a} = 'Lunge-Associated';
    else
        useableflukebeatdeep.MaxOrNormal{a} = 'Routine';
    end
    % Advance lunge of interest
    while useableflukebeatdeep.BeatStartFBD(a) > thislungeI && thislunge < numlunges
        thislunge = thislunge + 1;
        thislungeI = LungeI(thislunge);
    end
end
clear a;

for a = 1:length(useableflukebeatdeep.BeatStartFBD);
    if isnan(LungeI(1:end));
        useableflukebeatdeep.MaxOrNormal{a} = 'Unknown';
    end
end


%% %% Speed Difference Analysis - In response to reviewer comments on Gough et al. 2021 (Can Also use through MultiSelectHayden.m)
for a=1:length(useableflukebeatdeep.BeatStart);
    useableflukebeatdeep.BodyMass{a} = mass;
    useableflukebeatdeep.BeatLength{a} = (useableflukebeatdeep.BeatEnd(a)-useableflukebeatdeep.BeatStart(a))/10;
    useableflukebeatdeep.BeatStartSpd{a} = jigspeed(useableflukebeatdeep.BeatStart(a));
    useableflukebeatdeep.BeatEndSpd{a} = jigspeed(useableflukebeatdeep.BeatEnd(a));
    useableflukebeatdeep.SpdChange{a} = useableflukebeatdeep.BeatEndSpd{a}-useableflukebeatdeep.BeatStartSpd{a};
    useableflukebeatdeep.SpdChngPerc{a} = abs(useableflukebeatdeep.SpdChange{a})/useableflukebeatdeep.BeatStartSpd{a}*100;
    useableflukebeatdeep.DistTravel{a} = (0.5*(useableflukebeatdeep.SpdChange{a}/useableflukebeatdeep.BeatLength{a})*(useableflukebeatdeep.BeatLength{a}^2))+(useableflukebeatdeep.BeatStartSpd{a}*useableflukebeatdeep.BeatLength{a}); % equation given by Jean
    useableflukebeatdeep.ThrustDragOffset{a} = ((useableflukebeatdeep.BodyMass{a}/2)*((useableflukebeatdeep.BeatEndSpd{a}^2)-(useableflukebeatdeep.BeatStartSpd{a}^2)))/useableflukebeatdeep.DistTravel{a}; % equation given by Jean
end
clear a;


%% Calculating Thrust Power (Watts), Drag Coefficient, and Reynolds Number For Each Flukebeat (IF ANIMAL ALREADY HAS TAILBEATS IN Droned Tailbeats Info Hayden, ONLY RUN THIS SECTION)
yatesequations = xlsread('D:\REU Hayden\Digitizing Yates Figure\Digitized Yates Figure.xlsx',2,'A2:E82');
yatesequations2 = xlsread('D:\REU Hayden\Digitizing Yates Figure\Digitized Yates Figure.xlsx',4,'A2:E82');

density = 1000; % density of seawater
alpha = 0.523; % approximated as 30 degrees or 0.523 radians
h = 0.5*(0.2*TL); % heave amplitude (midline to peak) approximated as one-half of 20 percent of body length
Ct = 0;
efficiency = 0;

if strcmp(Species, 'Humpback') == 1;
    k = 0.05;
else
    k = 0.03;
end

for a=1:length(useableflukebeatdeep.BeatStart); % calculate out the coefficient of thrust (Ct) and efficiency using the graphs in Yates
    w = 2*pi*(useableflukebeatdeep.OsFreq(a));
    feathering = (((useableflukebeatdeep.AvgSpeeds(a))*alpha)/(w*h));
    round2 = @(x, d) round(x * 10^d)/10^d; 
    roundedfeathering = round2(feathering,2);
    reduced = ((w*C)/(useableflukebeatdeep.AvgSpeeds(a)));
    for b=1:81;
        roundedyates = round2((yatesequations(b,1)),2);
        if roundedfeathering == roundedyates;
            Ct = (((yatesequations(b,2))*(reduced^3)) + ((yatesequations(b,3))*(reduced^2)) + ((yatesequations(b,4))*reduced) + (yatesequations(b,5)));
            efficiency = (((yatesequations2(b,2))*(reduced^3)) + ((yatesequations2(b,3))*(reduced^2)) + ((yatesequations2(b,4))*reduced) + (yatesequations2(b,5)));
            break
        else if roundedfeathering >= 0.8000;
                Ct = (((yatesequations(81,2))*(reduced^3)) + ((yatesequations(81,3))*(reduced^2)) + ((yatesequations(81,4))*reduced) + (yatesequations(81,5)));
                efficiency = (((yatesequations2(81,2))*(reduced^3)) + ((yatesequations2(81,3))*(reduced^2)) + ((yatesequations2(81,4))*reduced) + (yatesequations2(81,5)));
            break
            else
                continue
            end
        end
    end
    useableflukebeatdeep.Eff(a) = efficiency;
    useableflukebeatdeep.Thrust(a) = (0.5*density*Ct*((useableflukebeatdeep.AvgSpeeds(a))^3)*Fa*((h/C)^2))*efficiency; % calculates the thrust power for each tailbeat
    useableflukebeatdeep.DragCoeffEqual(a) = ((useableflukebeatdeep.Thrust(a))/(0.5*density*Sa*((useableflukebeatdeep.AvgSpeeds(a))^3))); % calculates the drag coefficient for each tailbeat
    useableflukebeatdeep.DragCoeffReal(a) = ((useableflukebeatdeep.Thrust(a)-(useableflukebeatdeep.AvgSpeeds(a)*k*mass*((useableflukebeatdeep.BeatEndSpd{a}-useableflukebeatdeep.BeatStartSpd{a})/useableflukebeatdeep.BeatLength{a})))/(0.5*density*Sa*((useableflukebeatdeep.AvgSpeeds(a))^3)));
    useableflukebeatdeep.Reynolds(a) = ((TL*(useableflukebeatdeep.AvgSpeeds(a)))/(1.044*10^-6)); % calculates the reynolds number for each tailbeat
    useableflukebeatdeep.PercentOfThrust{a} = (useableflukebeatdeep.ThrustDragOffset{a}/useableflukebeatdeep.Thrust(a))*100;
end
clear a;

avgthrust = mean(useableflukebeatdeep.Thrust);
avgdrag = mean(useableflukebeatdeep.DragCoeffReal);
avgreynolds = mean(useableflukebeatdeep.Reynolds);

disp('Average Of Thrust Power For Each Flukebeat: ');
disp(avgthrust);
disp('Average Of Drag Coefficient For Each Flukebeat: ');
disp(avgdrag);
disp('Average Reynolds Number During Fluking: ');
disp(avgreynolds);

for a=1:length(useableflukebeatdeep.BeatStart); % adds in extra metadata to the spreadsheet that will be output (useableflukebeatdeep)
    useableflukebeatdeep.TotLength(a) = TL;
    useableflukebeatdeep.WettedSurfArea(a) = Sa;
    useableflukebeatdeep.FlukeArea(a) = Fa;
    useableflukebeatdeep.Species(a) = Species;
    useableflukebeatdeep.DeployID(a) = cellstr(INFO.whaleName);
end
clear a;

save(['C:\Users\William Gough\Documents\Academic Materials\Stanford University\Github Repositories\Hayden_REU_Project\Hayden_REU_Project-master\output data\' INFO.whaleName 'useableflukebeatdeep.mat'], 'useableflukebeatdeep');
writetable(useableflukebeatdeep,['C:\Users\William Gough\Documents\Academic Materials\Stanford University\Github Repositories\Hayden_REU_Project\Hayden_REU_Project-master\output data\' INFO.whaleName 'useableflukebeatdeep.csv']);

outfile = ['C:\Users\William Gough\Documents\Academic Materials\Stanford University\Github Repositories\Hayden_REU_Project\Hayden_REU_Project-master\output data\' INFO.whaleName 'useableflukebeatdeep.csv'];
if exist(outfile, 'file') == 2;
    delete(outfile);
end
writetable(useableflukebeatdeep,outfile);

