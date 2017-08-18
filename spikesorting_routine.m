% spike sorting routine for open-ephys data

addpath(genpath('C:\Users\Jake\Repos\KiloSort')) % path to kilosort folder
addpath(genpath('C:\Users\Jake\Repos\npy-matlab')) % path to npy-matlab scripts
addpath(genpath('C:\Users\Jake\Dropbox\MatlabCode\Repos\RigBuildPhotodiodeTest\analysis-tools-master'))
%%
fpath = uigetdir();


%%

makeAtlasChannelMap(fpath);

ops=buildConfigFile(fpath);
ops.artifactThresh = 10000;



if ops.GPU
    fprintf('Using GPU for faster processing\n')
    gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
end

%% artifact detection goes here?

% additional preprocessing steps before running kilosort
% 1) save high-pass filtered version of the binary file for visualization
%    in Phy
% 2) detect large events that span multiple channels and remove them
ops=removeArtifacts(ops);

%% check that the timing lines up
checkRawBinaryMatch(ops, 10e3)

%% main spike-sorting routine
[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)


%% saving
% save matlab results file
fprintf('saving matlab results file\n')
save(fullfile(ops.root,  'rez.mat'), 'rez', 'ops', '-v7.3');

rez                = merge_posthoc2(rez);
fprintf('saving python files for Phy\n')
% save python results file for Phy
rezToPhy(rez, ops.root);

% remove temporary file
delete(ops.fproc);