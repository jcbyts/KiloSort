function ops=buildConfigFile(fpath)


[~,thisComputer]=system('hostname');
switch thisComputer(1:end-1)
    case 'Gravedigger'
        ops.GPU         = 1; % whether to run this code on an Nvidia GPU (much faster, mexGPUall first)
        ops.parfor      = 1; % whether to use parfor to accelerate some parts of the algorithm
        ops.verbose     = 1; % whether to print command line progress		
        ops.showfigures = 0; % whether to plot figures during optimization
    otherwise
        ops.GPU         = 0; % whether to run this code on an Nvidia GPU (much faster, mexGPUall first)
        ops.parfor      = 0; % whether to use parfor to accelerate some parts of the algorithm
        ops.verbose     = 1; % whether to print command line progress		
        ops.showfigures = 1; % whether to plot figures during optimization
end

ops.root                = fpath;
ops.fbinary             = fullfile(ops.root, 'rawbinary.dat'); % will be created for 'openEphys'

% --- convert open ephys to binary

% check for multiple blocks
fl=dir(fullfile(ops.root, '*CH*.continuous'));
fl_=dir(fullfile(ops.root, '*CH*_*.continuous'));
if numel(fl_) > numel(fl)
    fl=fl_;
    str='CH(\w*)_';
else
    str='CH(\w*).';
end
    
ch=nan(numel(fl),1);
for j=1:numel(fl)
    str='CH(\w*).';
    ch(j)=cellfun(@(x) str2double(x{1}), regexp(fl(j).name,str, 'tokens'));
    if isnan(ch(j))
        str='CH(\w*)_';
        ch(j)=cellfun(@(x) str2double(x{1}), regexp(fl(j).name,str, 'tokens'));
    end
end

ops.NchanTOT=numel(unique(ch));
ops.Nchan=ops.NchanTOT;
for j = 1:numel(fl)
    fid=fopen(fullfile(ops.root,fl(j).name));
    fseek(fid,0,'bof');
    hdr=fread(fid,1024,'char*1');
    eval(char(hdr'))
    
    row = mod(j, ops.Nchan);
    if row==0
        row=ops.Nchan;
    end
    col = ceil(j/ops.Nchan);
    ops.chHeaders(row,col) = header;
end

sampleRates = arrayfun(@(x) x.sampleRate, ops.chHeaders(:));
sampleRate = unique(sampleRates);
assert(numel(sampleRate)==1, 'Different sampling rates found')
ops.fs=sampleRate;

if ~exist(ops.fbinary, 'file')
    ops = convertOpenEphysToRawBInary(ops);
end

ops.datatype            = 'dat';  % binary ('dat', 'bin') or 'openEphys'
ops.fproc               = fullfile(fpath, 'tempWh.dat'); % residual from RAM of preprocessed data				
% define the channel map as a filename (string) or simply an array
ops.chanMap             = fullfile(fpath, 'chanMap.mat'); % make this file using createChannelMapFile.m
if ~exist(ops.chanMap, 'file')
	ops.chanMap = 1:ops.Nchan; % treated as linear probe if unavailable chanMap file
end

ops.Nfilt               = 2*(ops.Nchan+32-mod(ops.Nchan,32));  % number of clusters to use (2-4 times more than Nchan, should be a multiple of 32)     		
ops.nNeighPC            = 12; % visualization only (Phy): number of channnels to mask the PCs, leave empty to skip (12)		
ops.nNeigh              = 16; % visualization only (Phy): number of neighboring templates to retain projections of (16)		
		
% options for channel whitening		
ops.whitening           = 'noSpikes'; % type of whitening (default 'full', for 'noSpikes' set options for spike detection below)		
ops.nSkipCov            = 1; % compute whitening matrix from every N-th batch (1)		
ops.whiteningRange      = 32; % how many channels to whiten together (Inf for whole probe whitening, should be fine if Nchan<=32)		
		
ops.criterionNoiseChannels = 0.2; % fraction of "noise" templates allowed to span all channel groups (see createChannelMapFile for more info). 		

% other options for controlling the model and optimization		
ops.Nrank               = 3;    % matrix rank of spike template model (3)		
ops.nfullpasses         = 6;    % number of complete passes through data during optimization (6)		
ops.maxFR               = 20000;  % maximum number of spikes to extract per batch (20000)		
ops.NotchFilter60       = false;
ops.softwareReferencing = 'none';
% ops.fshigh              = 200;   % frequency for high pass filtering		
% ops.fslow             = 2000;   % frequency for low pass filtering (optional)
ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection		
ops.scaleproc           = 200;   % int16 scaling of whitened data		
ops.NT                  = 128*1024+ ops.ntbuff;% this is the batch size (try decreasing if out of memory) 		
% for GPU should be multiple of 32 + ntbuff		
		
% the following options can improve/deteriorate results. 		
% when multiple values are provided for an option, the first two are beginning and ending anneal values, 		
% the third is the value used in the final pass. 		
% ops.Th               = [6 12 12];    % threshold for detecting spikes on template-filtered data ([6 12 12])		
% ops.lam              = [10 30 30];   % large means amplitudes are forced around the mean ([10 30 30])		
% ops.nannealpasses    = 4;            % should be less than nfullpasses (4)		
% ops.momentum         = 1./[20 1000];  % start with high momentum and anneal (1./[20 1000])		
% ops.shuffle_clusters = 1;            % allow merges and splits during optimization (1)		
% ops.mergeT           = .1;           % upper threshold for merging (.1)		
% ops.splitT           = .1;           % lower threshold for splitting (.1)

ops.Th               = [4 10 10];    % threshold for detecting spikes on template-filtered data ([6 12 12])		
ops.lam              = [5 5 5];      % large means amplitudes are forced around the mean ([10 30 30])		
ops.nannealpasses    = 4;            % should be less than nfullpasses (4)		
ops.momentum         = 1./[20 400];  % start with high momentum and anneal (1./[20 1000])		
ops.shuffle_clusters = 1;            % allow merges and splits during optimization (1)		
ops.mergeT           = .1;           % upper threshold for merging (.1)		
ops.splitT           = .1;           % lower threshold for splitting (.1)	
		
% options for initializing spikes from data		
ops.initialize      = 'fromData';    %'fromData' or 'no'		
ops.spkTh           = -6;      % spike threshold in standard deviations (4)		
ops.loc_range       = [3  1];  % ranges to detect peaks; plus/minus in time and channel ([3 1])		
ops.long_range      = [30  6]; % ranges to detect isolated peaks ([30 6])		
ops.maskMaxChannels = 5;       % how many channels to mask up/down ([5])		
ops.crit            = .65;     % upper criterion for discarding spike repeates (0.65)		
ops.nFiltMax        = 20000;   % maximum "unique" spikes to consider (10000)		
		
% load predefined principal components (visualization only (Phy): used for features)		
dd                  = load('PCspikes2.mat'); % you might want to recompute this from your own data		
ops.wPCA            = dd.Wi(:,1:7);   % PCs 		
		
% options for posthoc merges (under construction)		
ops.fracse  = 0.1; % binning step along discriminant axis for posthoc merges (in units of sd)		
ops.epu     = Inf;		

if ispc
    dmem=memory;
    ops.ForceMaxRAMforDat   = dmem.MaxPossibleArrayBytes;
else
    ops.ForceMaxRAMforDat   = 20e9; % maxim
end