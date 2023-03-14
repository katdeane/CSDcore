% this is the most basic csd generation script which can then be inserted
% into a loop or pipeline for generating multiple CSDs. Have a look at
% SingleTrialCSD* scripts in my pipelines for a fuller version and for more
% types of output. Email me at katrinad@ucr.edu for assistance.

clear; % clear your workspace
% make sure Matlab is open in the main current directory
% cd('E:\CSDcore') % change directory to the main directory in your folders
homedir = pwd; % print working directory
addpath(genpath(homedir)); % genpath includes all subfolders into working directory now

% load or generated LFP data (channels x time x trials)
LFP = rand(16,1000,20); % generate a 16 channel probe over 1000 ms for 20 trials
%load('LFP_data.mat','LFP'); % load the mat with LFP data 

% assuming (channel x time x trials) data:
Ntrials = size(LFP,3); % determine the number of trials 

%% CSD
% Parameters for filtering
chan_dist = 50;    % channel distance in your electrodes in µmh
kernel    = 450;   % 600 is very smooth, 300 is more accurate, 450 is happy medium
kernelxchannel = kernel./chan_dist; % kernel size in terms of number of channels
hammsiz   = kernelxchannel+(rem(kernelxchannel,2)-1)*-1; % linear extrapolation and running avg
paddsiz   = floor(hammsiz/2)+1; % padding size necessary 

% preallocate:
singletrialCSD = NaN(size(LFP));
% loop through trials
for itrial = 1:Ntrials
    curLFP = LFP(:, :, itrial);                            % current LFP trial
    paddit = padd_linex(curLFP,paddsiz,1);                 % create padding
    hammit = filt_Hamm(paddit,hammsiz,1);                  % apply hamming filter
    singletrialCSD(:,:,itrial) = (get_csd(hammit,1,chan_dist,1))*10^3;  % make csd
end

% you can take the average CSD for faster calculation if you don't need
% single trial 
CSD = nanmean(singletrialCSD,3);

figure; % plot CSD
imagesc(CSD) % this command displays image with scaled colors
caxis([-0.001 0.001]) % change the color axis - play with this based on your scaling, see Deane et al 2020 for example CSD output visual
colormap jet % change the colormap to jet for CSDs (perula is default)
colorbar % include the colorbar (and scale) into the image

%% Average rectified CSD 
% rectifying just means taking the absolute values 
% This gives you the overall activity in the column over time 

singletrialAVREC = nanmean(abs(singletrialCSD)); % computes rectified mean (see Schroeder et al. Cereb Cortex epub ahead of print September 2006)
AVREC = nanmean(abs(CSD));

figure; % plot AVREC trace
plot(AVREC)



%% Relative residuals CSD
% the idea here is to see roughly what percentage must be coming from the
% thalamus (vertical input causing balanced sinks around the electrode) and
% what must be from cortial sources (horizontal input causing unbalanced
% sinks around the electrode). 

% preallocate:
singletrialRELRES = NaN(Ntrials,size(CSD,2));
for itrial = 1:Ntrials
    curCSD = singletrialCSD(:, :, itrial);
    % compute absolute residues 
    sumCSD = sum(curCSD);           
    % compute relative residues
    singletrialRELRES(itrial,:) = sumCSD./sum(abs(curCSD));  
end

RELRES = nanmean(singletrialRELRES,1);

figure; % plot RELRES trace
plot(RELRES)

%% 

% Analysis can be further broken down into layer trace analyses and sink activity analyses. 
% Above are the broadest strokes for CSD for cortical data from laminar probes.

