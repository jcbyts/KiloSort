function ops = convertOpenEphysToRawBInary(ops)

path_to=fileparts(ops.fbinary);
if isempty(path_to)
    fname       = fullfile(ops.root, sprintf('%s.dat', ops.fbinary));
else
    fname       = ops.fbinary;
end
fidout      = fopen(fname, 'w');
%
clear fs
fl=dir(fullfile(ops.root, '*CH*.continuous'));
fl_=dir(fullfile(ops.root, '*CH*_*.continuous'));
if numel(fl_) > numel(fl)
    for j = 1:ops.Nchan
            fs{j} = dir(fullfile(ops.root, sprintf('*CH%d_*.continuous', j) )); % if separate files are saved by open ephys gui
    end
else
    for j = 1:ops.Nchan
        fs{j} = dir(fullfile(ops.root, sprintf('*CH%d.continuous',j) )); % if single files are saved
    end
end

nblocks = cellfun(@(x) numel(x), fs);
if numel(unique(nblocks))>1
   error('different number of blocks for different channels!') 
end
%
nBlocks     = unique(nblocks);
nSamples    = 1024;  % fixed to 1024 for now!

fid = cell(ops.Nchan, 1);

tic
for k = 1:nBlocks
    for j = 1:ops.Nchan
        fid{j}             = fopen(fullfile(ops.root, fs{j}(k).name));
        % discard header information
        fseek(fid{j}, 1024, 0);
    end
    %
    nsamps = 0;
    flag = 1;
    while 1
        samples = zeros(nSamples * 1000, ops.Nchan, 'int16');
        for j = 1:ops.Nchan
            collectSamps    = zeros(nSamples * 1000, 1, 'int16');
            
            rawData         = fread(fid{j}, 1000 * (nSamples + 6), '1030*int16', 10, 'b');

            nbatches        = ceil(numel(rawData)/(nSamples+6));
            for s = 1:nbatches
                rawSamps = rawData((s-1) * (nSamples + 6) +6+ [1:nSamples]);
                collectSamps((s-1)*nSamples + [1:nSamples]) = rawSamps;
            end
            samples(:,j)         = collectSamps;
        end
        
        if nbatches<1000
            flag = 0;
        end
        if flag==0
            samples = samples(1:s*nSamples, :);
        end
       
        samples         = samples';
        fwrite(fidout, samples, 'int16');
        
        nsamps = nsamps + size(samples,2);
        
        if flag==0
            break;
        end
    end
    ops.nSamplesBlocks(k) = nsamps;
    
    for j = 1:ops.Nchan
       fclose(fid{j}); 
    end
    
end
    
fclose(fidout);

toc