function makeAtlasChannelMap(fpath)
% create a channel Map file for 32 Channel Atlas probe

% here I know a priori what order my channels are in.  So I just manually 
% make a list of channel indices (and give
% an index to dead channels too). chanMap(1) is the row in the raw binary
% file for the first channel.



% the first thing Kilosort does is reorder the data with data = data(chanMap, :).

% 32Chan Atlas Acute with drivemoutn
chanMap = [5 21 7 23 3 19 12 28 14 30 8 24 6 22 10 26 9 25 4 20 11 27 16 32 13 29 2 18 1 17 15 31];

% 32chan Atlas Chronic ZIF to omnetics
% chanMap = [9    26     8    24    10    25     6    23    11    27     7    22    12    28     4    21    13    29, ...
%     5    20    14    30     2    19    15    31     3    18    32    17    16     1];
chanMap = [8 24 9 25 7 23 10 26 11 27 12 28 5 21 4 20 3 19 2 18 1 6 13 32 14 31 15 30 16 29 17 22];

% El = 1:32;
% ZIF = [18 38 23 3 17 37 24 4 16 36 25 5 15 35 26 6 14 34 27 7 13 33 28 8 12 32 29 9 30 10 11 31];
% 
% Omn = [nan nan  35:-1:20 nan nan nan nan 2:17 nan];
% intan = [nan 8 6 7 4 5 2 3 32 1 31 30 29 28 27 25 26 nan, ...
%     nan 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 nan];



% Now we declare which channels are "connected" in this normal ordering, 
% meaning not dead or used for non-ephys data
connected = true(32, 1); % connected(1:2) = 0;

% now we define the horizontal (x) and vertical (y) coordinates of these
% 34 channels. For dead or nonephys channels the values won't matter. Again
% I will take this information from the specifications of the probe. These
% are in um here, but the absolute scaling doesn't really matter in the
% algorithm. 

xcoords = zeros(1,32);
ycoords = 50*(1:32); 

% Often, multi-shank probes or tetrodes will be organized into groups of
% channels that cannot possibly share spikes with the rest of the probe. This helps
% the algorithm discard noisy templates shared across groups. In
% this case, we set kcoords to indicate which group the channel belongs to.
% In our case all channels are on the same shank in a single group so we
% assign them all to group 1. 

kcoords = ones(1,32);

% at this point in Kilosort we do data = data(connected, :), ycoords =
% ycoords(connected), xcoords = xcoords(connected) and kcoords =
% kcoords(connected) and no more channel map information is needed (in particular
% no "adjacency graphs" like in KlustaKwik). 
% Now we can save our channel map for the eMouse. 

% would be good to also save the sampling frequency here
fs = 32000; 

save(fullfile(fpath, 'chanMap.mat'), 'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords', 'fs')