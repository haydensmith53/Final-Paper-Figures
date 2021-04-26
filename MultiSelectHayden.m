%% Run Tailbeat Data From All Individuals (Or A Chosen Set)
%  Imports tailbeat data from each individual and runs Chapter 1 scripts.
clear all;

[filename, filelocs] = uigetfile('MultiSelect','on');
for a=1:length(filename);
    load([filelocs filename{a}]);
    FlukingKinoHydroFinalSmith2021Update;
    clearvars -except filename filelocs morphometrics
end
clear a;

%% Combine All Data Together
[filename, filelocs] = uigetfile('*useableflukebeatdeep.csv', 'MultiSelect','on');
AllDronedFlukebeats = readtable([filelocs filename{1}]); 
for a=2:length(filename);
    thisWhale = readtable([filelocs filename{a}]);
    AllDronedFlukebeats = vertcat(AllDronedFlukebeats, thisWhale);
end
clear a;

outfile = ['C:\Users\William Gough\Documents\Academic Materials\Stanford University\Github Repositories\Final-Paper-Figures\AllDronedFlukebeatsFinalized.csv'];
if exist(outfile, 'file') == 2;
    delete(outfile);
end
writetable(AllDronedFlukebeats,outfile);

%% Perform Speed Change Analysis (Not Necessary At This Stage)
AllData = readtable('C:\Users\William Gough\Documents\Academic Materials\Stanford University\Github Repositories\Final-Paper-Figures\AllDronedFlukebeatsFinalized.csv'); % read in all of the morphometric data
morphometrics = readtable('C:\Users\William Gough\Documents\Academic Materials\Stanford University\Github Repositories\Final-Paper-Figures\Finalized Data Sheet For Hayden.xlsx'); % read in all of the morphometric data

[filename, filelocs] = uigetfile('MultiSelect','on');
for a=1:length(filename);
    load([filelocs filename{a}]);
    SpeedChangeFinalSmith2020;
    clearvars -except filename filelocs AllData morphometrics
end
clear a;

outfile = ['C:\Users\William Gough\Documents\Academic Materials\Stanford University\Github Repositories\Final-Paper-Figures\AllDronedFlukebeatsFinalized.csv'];
if exist(outfile, 'file') == 2;
    delete(outfile);
end
writetable(AllData,outfile);