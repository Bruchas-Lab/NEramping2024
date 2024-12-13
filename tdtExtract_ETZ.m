 % Code written by Christian Pedersen
% Michael Bruchas Lab - UW

% Extracts photometry data from TDT Data Tank and outputs time series as
% .mat file

% REQUIRES tdt2mat.m file to be in same folder (filepath)

%% Reset MatLab workspace - clears all variables and the command window

clear all;  % clear all variables
close all;  % close all open graphs
clc   % clear command window


%%

%%%%%%%%%%%%%%%%%%%%%%  EDIT THESE FIELDS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% in path (white bar above) navigate to D:\MATLAB\Data and make new folder for condition   

% point to tanks (enter file path)
path_to_data = 'D:\Photometry\CFDShortDay3B\Tanks';

% set up Tank name, variables to extract
tankdir = [path_to_data];
% tankname = 'Eric_stim_temp-230824-135010'; % name of your tank
% tankname = 'Avi_stim-230330-145256';
% tankname = 'Avi_stim-230801-132002';
tankname = 'Avi_stim-240105-123944';
blockname = '392-4_Day3B-240816-131050'; % name of your file

% pick any file name to save time series (must end in .mat)
filename = blockname;
%after changing path_to_data, tankname, blockname, run code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% %Bruchas cart:
storenames = {'470A'}; % name of stores to extract from TDT (usu. 4-letter code) 
% %LMag is the demodulated data, may also have other timestamps etc
% 
storenames2 = {'405A'};

% Static system Box C:
% storenames = {'465C'}; % name of stores to extract from TDT (usu. 4-letter code) 
%LMag is the demodulated data, may also have other timestamps etc
%
% storenames2 = {'405C'};



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract

for k = 1:numel(storenames)
  storename = storenames{k};
  S{k} = tdt2mat(tankdir, tankname, blockname, storename);
end

for k = 1:numel(storenames2)
  storename2 = storenames2{k};
  S2{k} = tdt2mat(tankdir, tankname, blockname, storename2);
end

% Massage data and get time stamps

LMag = S{1}; %add more if you extracted more stores above
% LMag2 = S{2};
% For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
chani1 = LMag.channels==1;
chani2 = LMag.channels==2;

LMag2 = S2{1}; %add more if you extracted more stores above
% LMag2 = S{2};
% For 2-color rig, LMag data is on channels 1 and 2, channel 1 = 470nm, channel 2 = 405nm
chani21 = LMag2.channels==1;
chani22 = LMag2.channels==2;
% chani21 = LMag2.channels==1;
% chani22 = LMag2.channels==2;

% Get LMag data as a vector (repeat for each channel)
rawdat1 = LMag.data(chani1,:);
rawdat1 = reshape(rawdat1', [],1); % unwrap data from m x 256 array
% dat2 = LMag.data(chani21,:);
% dat2 = reshape(dat2', [],1); % unwrap data from m x 256 array

% Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
ts = LMag.timestamps(chani1);
t_rec_start = ts(1);

ts = ts-ts(1); % convert from Unix time to 'seconds from block start'
ts = bsxfun(@plus, ts(:), (0:LMag.npoints-1)*(1./LMag.sampling_rate));
ts = reshape(ts',[],1);

%%%%%%%%%%%%%%%%%%


dat2 = LMag2.data(chani21,:);
dat2 = reshape(dat2', [],1); % unwrap data from m x 256 array
% dat2 = LMag.data(chani21,:);
% dat2 = reshape(dat2', [],1); % unwrap data from m x 256 array

% Get LMag timestamps (use chani1 - timestamps should be the same for all Wpht channels
ts2 = LMag2.timestamps(chani21);
t_rec_start2 = ts2(1);

ts2 = ts2-ts2(1); % convert from Unix time to 'seconds from block start'
ts2 = bsxfun(@plus, ts2(:), (0:LMag2.npoints-1)*(1./LMag2.sampling_rate));
ts2 = reshape(ts2',[],1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subtracted 405 signal from 470 signal % raw final signal
% rawdat1 = rawdat1(1:202368);
% ts = ts(1:202368);
subdat = rawdat1-dat2;
Dts=ts;
% plot raw 405 ch, raw 470 ch and subtracted signal
figure(2)
hold on
plot(ts,rawdat1,'b');
plot(ts2,dat2,'r');
title('both signals')
plot(ts,subdat,'g');
xlabel('time(s)')
ylabel('amplitude')
% axis([0 3600 160 210])
hold off


%% Bandpass filter data 
Fs = round(1/(max(Dts)/length(Dts))); %sample frequency Hz

%% IF PHOTOMETRY RECORDING BASELINE DROPS OFF:

%7200 cutoff P6-3 PR 1
% subdat = subdat(1:Fs*7200);
% Dts = Dts(1:length(subdat));

% % 2800 cutoff P11-1 PR 2
% subdat = subdat(1:Fs*2850);
% Dts = Dts(1:length(subdat));

%% instead of HPF: fit 2nd order exponential and subtract

%f2 = fit(Dts,subdat,'exp2');

%fitcurve= f2(Dts);


p1 = polyfit(Dts,subdat,4);
p2 = polyfit(Dts,rawdat1,4);
%Evaluate the polynomial on a finer grid and plot the results.

fitcurve = polyval(p1,Dts);
fitcurve2 = polyval(p2,Dts);


dataFilt = subdat - fitcurve;

dataFilt2 = rawdat1 - fitcurve2;

%dataFilt = subdat;

% f2 = fit(Dts,dat2,'exp2');
% fitcurve = f2(Dts);
% dataFilt = rawdat1 - fitcurve;

%% Calculate dF/F: this is my method,isosbestos and 465 calculated independently to derive dF v
%%%Subtract normalized raw 470 signal
        GEf1 = fit(Dts,rawdat1,'poly2');
        GEfitcurve1 = GEf1(Dts);
        GE470fit = rawdat1-GEfitcurve1;
        dF470 = (GE470fit./GEfitcurve1);
    %%%Subtract normalized raw 405 signal
        GEf2 = fit(Dts, dat2,'poly2');
        GEfitcurve2 = GEf2(Dts);
        GE405fit = dat2-GEfitcurve2;
        dF405 = (GE405fit./GEfitcurve2);
        
%%%Calculate dF/F
        GEF = dF470-dF405;  %%%Subtract normalized 405 from 470 signal (dF470-dF405) for F 
        GEf3 = fit(Dts,GEF,'poly2'); %%%fit curve to F
        GEF0 = GEf3(Dts);    %%F0=F(fitted)
        GEdata1 = (GEF-GEF0)/abs(median(dF470)); %%(dF470-dF405)/absolute value of (median(dF470))


%% low pass filter

% N = 1;  % order
% cutoff_Hz = 10;  %
% [b,a]=butter(N,cutoff_Hz/(Fs/2),'low');
% 
% dataFilt = filter(b,a,datatemp1);

%% calculate dF/F

normDat1 = (dataFilt - median(dataFilt))./abs(median(rawdat1)); % this gives deltaF/F
normDat = normDat1.*100; % make a percentage

normDat2 = (dataFilt2 - median(dataFilt2))./abs(median(rawdat1));
normDat2 = normDat2.*100;
data1 = normDat;
data2 = normDat2;

% plot subtracted signal w/o baseline correction
figure(6);
plot(Dts,subdat,Dts,fitcurve);
title('Subtracted signal')
xlabel('time(s)')
ylabel('raw F')

% plot final, baseline-corrected signal
figure(7);
plot(Dts,data1);
title('Signal with corrected baseline')
xlabel('time(s)')
ylabel('deltaF/F')
% ylim([-5 5])

% plot final, baseline-corrected signal
figure(8);
plot(Dts,data2);
title('Signal with corrected baseline')
xlabel('time(s)')
ylabel('deltaF/F')
% ylim([-20 20])

%% plot GEdata1 (different method of df/f)
figure(9);
plot(Dts,GEdata1);



%% Save file as .mat file with specified filename

save(filename, 'Dts','rawdat1','subdat','data1','data2');











