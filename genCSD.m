
homedir = pwd;
addpath(genpath(homedir));
        
% this is the most basic csd generation script which can then be inserted
% into a loop or pipeline for generating multiple CSDs. Have a look at
% SingleTrialCSD* scripts in my pipelines for a fuller version and for more
% types of output. 

% load or generated LFP data (channels x time x trials)
LFP = rand(16,1000,20); % a 16 channel probe over 1000 ms for 20 trials

Ntrials = size(LFP,3); % variable for number of trials

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
    curLFP = LFP(:, :, itrial);                           % current LFP trial
    paddit = padd_linex(curLFP,paddsiz,1);                 % create padding
    hammit = filt_Hamm(paddit,hammsiz,1);                  % apply hamming filter
    singletrialCSD(:,:,itrial) = (get_csd(hammit,1,chan_dist,1))*10^3;  % make csd
end

% you can take the average CSD for faster calculation if you don't need
% single trial 
CSD = mean(singletrialCSD,3);

% Average rectified CSD gives you the overall activity in the column over
% time 
singletrialAVREC = mean(abs(singletrialCSD)); % computes rectified mean (see Schroeder et al. Cereb Cortex epub ahead of print September 2006)
AVREC = mean(abs(CSD));

%% plotting

figure; % plot CSD
imagesc(CSD)
colormap('jet')
colorbar

figure; % plot AVREC trace
plot(AVREC)