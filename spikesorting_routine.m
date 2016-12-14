% spike sorting routine for open-ephys data
fpath = uigetdir();

makeAtlasChannelMap(fpath);

ops=buildConfigFile(fpath);

if ops.GPU     
    gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
end

%% artifact detection goes here?

% additional preprocessing steps before running kilosort
% 1) save high-pass filtered version of the binary file for visualization
%    in Phy
% 2) detect large events that span multiple channels and remove them
ops=removeArtifacts(ops);


%% main spike-sorting routine
[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)
rez                = merge_posthoc2(rez);

%% saving
% save matlab results file
fprintf('saving matlab results file\n')
save(fullfile(ops.root,  'rez.mat'), 'rez', '-v7.3');

fprintf('saving python files for Phy\n')
% save python results file for Phy
rezToPhy(rez, ops.root);

% remove temporary file
delete(ops.fproc);