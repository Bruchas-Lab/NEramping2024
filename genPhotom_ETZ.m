% general photometry processing script for preprocessed TDT files
% preprocessed .mat files should contain "Dts", "data1" 
% also requires txt file with times with event labels
% (1,2,...) separated by a space. One such file is included in the
% photometry package folder. The current code expects a similar format.
%
% To use: search CUSTOMIZE to find the lines in the script you need to adjust
%
% created by Christian Pederson 
% adjusted for use by Avi Matarasso, 02/21/21
% please email akmat@uw.edu with any questions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

codeDir = 'D:\MATLAB'; %CUSTOMIZE
dataDir = 'D:\MATLAB\Data\CFDDay1B'; %CUSTOMIZE

addpath(codeDir);

path_to_data =  {dataDir}; %CUSTOMIZE
% addpath('F:\code\PhotometryPackage');

    
%% CUSTOMIZE this whole section
%do you only want to look at the averaged data of *some* of the subjects?
% if you have subjects 0,1,2,4 in the directory you're working in, 
% and you want to look at 1 and 4, type [2 4];

%%%% commonStr == common string between your experiments of interest
region = 'grabNE'; 
commonStr = region; %regionON = 0; %mark as 0 if region not included in name
savename = 'CFD Day 1B';
titleName = savename;

rawtimings = dlmread('TIMING - CFD.txt'); %change to correct event times txt file
%%%%Txt files found at D:\MATLAB\PhotometryCode, for individual timings seeline 124

% timeBefore = 30;
timeBefore = 180;
% timeAfter  = 150;
timeAfter  = 17;
xlims = [-timeBefore timeAfter];
    
%%%%%%%% change the savename to get each timelocked trial (usually of stim)
doYouWantAll = 0; % make 0 if you want to only choose some of the subjects %and define which ones you want through whichSub
% if doing multiple paths, want to be careful and check whichSub you choose
whichSub = [1,3,4,5,7,8,9,10,11,13,14,15]; % which subjects do you want to look at/for multiple regions include all

%%% events what event do you want to look at
eventLabel  = 1; %i usually do 2 for tail lift and 1 for stim
manualEvent = 0; %1 if the scoring was done manually, 0 if TTL

%How many trials and do you want all?
allTrials = 1; %if you don't want all trials, make allTrials = 0
if ~allTrials 
    trials   = 3; % how many of eventLabel# did you do?
    trialsOI = 1:2; % 1:10 for example, which trials are you interested in
end

%%% do you want to z-score? if yes, make zQuest = 1, if no, make zQuest = 0
zQuest = 1;

trialNumber = 1; %adds one each iteration of session
normalQ = 0; %do you want to normalize to time immediately prior to event

%%% do you want each subject's timelocked and averaged traces and heat map?
% if averageQuest = 1, they will be shown when you run the code
averageQuest = 0;

%%% how much time after the event do you want to look at 
%timewindow = 30; % +/- seconds from event %usually 120
%specialCase = 2; %multiplier for timewindow to end. [-timewindow specialCase*timewindow]
% If only want to see time within +/- timewindow, keep specialCase as 1

%% Define Directory and constant variables

workdir = defineDir(path_to_data,commonStr);

%Define the working directories based on whichSub
if ~doYouWantAll
    workdir = workdir(whichSub);
end

% baselineTime = 180;
baselineTime = 15;
baselineSD = zeros(length(workdir),1);
baselineMu = zeros(length(workdir),1);

subName = cell(1, length(workdir));
allTimings    = cell(length(workdir),1);
FS = 1017.25; % sampling frequency of synapse
lengthOfBL = floor(abs(timeBefore)*FS);
lengthOfData = floor((timeAfter+timeBefore)*FS);

%CUSTOMIZE
dsFactor = 500; % how much do you want to downsample by?

proc = ['Data']; %add any processing steps
%ex: 'znData' would be z scored and normalized

%% Set up events and align them

for subjectIdx = 1:length(workdir) 

    clearvars -except sess k itchStack workdir savename trialname ...
         tempPath path_to_data sTs tT whichSub zQuest baselineTime baselineSD baselineMu...
        trialNumber subName subStr averageQuest firstQuest dsFactor trialsOI ...
        eventLabel allTrials timewindow specialCase normalQ lengthOfBL lengthOfData ...
        timeBefore timeAfter proc time subjectIdx manualEvent FS bl nData timingIdxs allData titleName xlims rawtimings
    
    %scroll through each photometry folder 
    cd(workdir(subjectIdx).folder)
    photoname = workdir(subjectIdx).name;
    subStr = [photoname(1:3) '-']; %CUSTOMIZE
    subNumb   = photoname(strfind(photoname, subStr) + length(subStr)); % will get sub
    load(photoname)
    
    %%% Make sure the Timings you use are consistent, or change them every time
    subjectLabel = photoname(1:5); %CUSTOMIZE based on your labels
    subName{subjectIdx} = subjectLabel;
%     txtName    = ['TIMING - ' subjectLabel '.txt'];
%     rawtimings = dlmread(txtName); %CUSTOMIZE
   
    
    eventTime = rawtimings(:,1);
    eventtype = round(rawtimings(:,2));
    eventTime = eventTime(eventtype==eventLabel);  
    if manualEvent
        tank = regexp(workdir(subjectIdx).folder,'\','split'); tank = tank{end};
        correctTime = max(Dts);
        vidDir = photoname(1:end-4); cDir = pwd;
        cd(vidDir); 
        vid = VideoReader([tank '_' vidDir '_Cam1.avi']); vidTime = vid.Duration;
        vidRatio = vidTime/correctTime;
        cd(cDir)
        allTimings{subjectIdx} = eventTime*vidRatio; %CUSTOMIZE
    else
        allTimings{subjectIdx} = eventTime;
    end
    
    %find baseline for each session
    baselineMu(subjectIdx) = mean(data1(1:ceil(baselineTime*FS)));
    baselineSD(subjectIdx) = std(data1(1:ceil(baselineTime*FS)));
    
    % Z score photometry data
    if zQuest
        data1 = (data1-mean(data1))./std(data1);
        proc = ['z' proc];
        warning('You z-scored your data')
    else
        warning('NOT z-scored')
    end    
    
    %Align the data to events 
    [nData, timingIdxs] = alignEvent(data1, FS, eventTime, timeBefore, timeAfter);

    proc = ['n' proc]; %update processing steps
    time = linspace(1/FS,length(data1)/FS,floor(length(data1)));
    numActualTrials  = size(nData,2);
    
    if allTrials
        trials   = numActualTrials; % how many of eventLabel# did you do?
        trialsOI = 1:trials; % 1:10 for example, which trials are you interested in
    end
    
    if numActualTrials < length(trialsOI) 
        cols = (subjectIdx-1)*length(trialsOI)+trialsOI;
        cols = cols(1:numActualTrials);
    else
        cols = (subjectIdx-1)*length(trialsOI)+trialsOI;
    end
    allData(:,cols) = nData(:,trialsOI); 
        
    %plot individual averages and heat map
    figureN = subjectIdx+100;
    plotIndividual(figureN, eventTime, data1, FS, subName, subjectIdx, dsFactor, time)
    
    if averageQuest
        heatFig = 200+subjectIdx;
        avgFig = 300+subjectIdx;
        plotHeatAndIndAvg(nData, heatFig, avgFig, timeBefore, timeAfter, subName, subjectIdx, dsFactor,FS)
    end
    
    if normalQ
        [cData,ncData] = centAndNormData(nData, FS, timeBefore); %normalized to pre-event
        NcData = cData/std(data1(1:round(180*FS))); 
    end
end

%% plot first stim session (not necessary)


%% save variables in folder
%{
if isfolder('varsAndFigs')
    cd('varsAndFigs')
    if exist(savename,'var')
        save(savename, 'allData','timingIdxs','sem','path_to_data','-append')
    else
        save(savename, 'allData','timingIdxs','sem','path_to_data')
    end
    cd('..')
else
    mkdir('varsAndFigs')
    cd('varsAndFigs')
    if exist(savename,'var')
        save(savename, 'allData','timingIdxs','sem','path_to_data','-append')
    else
        save(savename, 'allData','timingIdxs','sem','path_to_data','path_to_data')
    end
    cd('..')
    % close all
end
%}

%%
%shift everything based on baseline
%custom code for individual traces
for i = 1:size(allData,2)
    allDataDS(:,i) = decimate(allData(:,i),dsFactor);
end

% event triggered average
figure(61)

hold on
x = decimate(linspace(-timeBefore,timeAfter,length(nData)),dsFactor);
if any(any(isnan(allData)))
    cols2getridof = isnan(allData(1,:));
    allData(:,cols2getridof) = [];
    warning('some columns were deleted!')
end
y = mean(allData,2); 
centerY = mean(mean(allData(1:round(timeBefore*FS),:)));
y = (y - centerY); %center the data
y = decimate(y,dsFactor);
y = smoothdata(y,'SmoothingFactor',0.005);
sem = std(allData,0,2)./sqrt(size(allData,2)); % sem = std/sqrt(n)
eb = decimate(sem',dsFactor);
lineProps.col{1} = [0 0.5 0];

[M,I] = max(y);
adjust = mean(y(1:61));
M = M - adjust;
S = eb(I);
y = y - adjust;

mseb(x,y,eb,lineProps,1);
%plot(linspace(-timewindow,timewindow,length(photoPerLick)),photoPerLick,'b')%,linspace(-timewindow,timewindow,length(photoPerLick)),LickTrig,'y')
%legend('SEM','Mean')
L = line([0 0],[-3 20]);
%L = line([30 30],[-3 5]);
if exist('xlims','var')
    xlim(xlims)
else
    xlim([-20 timeAfter])
end
% ylim([-3 3])
ylim([-2 5])
%axis([-timewindow timewindow -0.6 1.2])
set(L,'Color','black')
xlabel('Time aligned to event (s)')
ylabel('z-score normalized to baseline')
title(titleName)
hold off
set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry

%% heat map (x dim: time, ydim: event, zdim: deltaF/F)
figure(62)

%uncomment these if you prefer to get 
%[~,idx1] = sort(mean(LickTrig1(:,90000:150000),2),'ascend');
%LickTrig1 = LickTrig1(idx1,:);

trials = length(trialsOI);
for i = 1:size(allData,2)
    adjust2 = mean(allData(1:30518,i));
    allData(:,i) = allData(:,i) - adjust2;
end
hold on
timevec=linspace(-timeBefore,timeAfter,length(allData));
imagesc(timevec,1:size(allData,2),allData')
L = line([0 0],[0 size(allData,1)+1]);
set(L,'Color','white')
set(L,'LineWidth',1.5)
% plot horizontal lines
for i = 1:length(workdir)
L2 = line([-150 size(allData,1)],[trials*i + 0.5 size(allData,2)+1]);
set(L2,'Color','white')
set(L2,'LineWidth',1.5)
end

xlabel('Time aligned to event (s)')
ylabel('Trial by mouse, bottom displays earliest trial')
%ytick(1:size(LickTrig,1))
title(titleName)
cb = colorbar;
title(cb,'Z score')
caxis([-2 2])
colormap(flipud(brewermap([],'YlGnBu')))
if exist('xlims','var')
    xlim(xlims)
else
    xlim([-20 100])
end
ylim([0.5 size(allData,2)+0.5])
set(gca, 'yticklabel', subName, 'ytick', trials/2 - 0.5 + (1:length(trialsOI):length(trialsOI)*length(subName)));
%for only one stim 
%set(gca, 'yticklabel', YLabelStrings, 'ytick', (1:length(YLabelStrings)));
hold off



%% Save All Figures
answer = questdlg('Would you like to save figures?', ...
	'Saving?', ...
	'Yes','No','No');
switch answer
    case 'Yes'
        figureSaveName = savename;
        if exist('varsAndFigs','dir')
            cd('varsAndFigs')
        else
            mkdir('varsAndFigs')
            cd('varsAndFigs')
        end
        tempdir = pwd;
        FolderName = tempdir;   % Your destination folder
        FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
        savefig(FigList, fullfile(FolderName,[figureSaveName '.fig']));
        save([savename '.mat']);
        disp('Figures have been saved!')        
        cd('..')
    case 'No'
        disp('You may manually save figures if you want.')
end

set(gcf, 'Renderer', 'Painters'); % for making sure the svg files don't come out blurry
