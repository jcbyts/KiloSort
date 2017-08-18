function ops = removeArtifacts(ops)
%REMOVEARTIFACTS high-pass filters and detects artifacts in
%electrophysiology traces
sizeThresh=getOr(ops, {'artifactThresh'}, 100); % artifacts must cross this threshold
chanThresh=getOr(ops, {'artifactNchans'}, 12); % artifact must span at least this number of channels

path_to=fileparts(ops.fbinary);
if isempty(path_to)
    fname       = fullfile(ops.root, sprintf('%s.dat', ops.fbinary));
else
    fname       = ops.fbinary;
end
[path_to, file, ext]=fileparts(ops.fbinary);
fnameOut=fullfile(path_to, [file '_hp' ext]);

fid      = fopen(fname);
fidout   = fopen(fnameOut, 'w');

NchanTOT = ops.NchanTOT;
NT       = ops.NT ;

d = dir(ops.fbinary);
ops.sampsToRead = floor(d.bytes/NchanTOT/2);

if ispc
    dmem         = memory;
    memfree      = dmem.MemAvailableAllArrays/8;
    memallocated = min(ops.ForceMaxRAMforDat, dmem.MemAvailableAllArrays) - memfree;
    memallocated = max(0, memallocated);
else
    memallocated = ops.ForceMaxRAMforDat;
end
nint16s      = memallocated/2;

NTbuff      = NT + 4*ops.ntbuff;
Nbatch      = ceil(d.bytes/2/NchanTOT /(NT-ops.ntbuff));
Nbatch_buff = floor(4/5 * nint16s/ops.Nchan /(NT-ops.ntbuff)); % factor of 4/5 for storing PCs of spikes
Nbatch_buff = min(Nbatch_buff, Nbatch);
ibatch = 0;
%% load data into patches, filter, compute covariance
if isfield(ops,'fslow')&&ops.fslow<ops.fs/2
    [b1, a1] = butter(3, [ops.fshigh/ops.fs,ops.fslow/ops.fs]*2, 'bandpass');
else
    [b1, a1] = butter(3, ops.fshigh/ops.fs*2, 'high');
end

isproc = zeros(Nbatch, 1);
while 1
    ibatch = ibatch + ops.nSkipCov;
    
    offset = max(0, 2*NchanTOT*((NT - ops.ntbuff) * (ibatch-1) - 2*ops.ntbuff));

    fseek(fid, offset, 'bof');
    buff = fread(fid, [NchanTOT NTbuff], '*int16');
    
    %         keyboard;
    
    if isempty(buff)
        break;
    end
    nsampcurr = size(buff,2);
    if nsampcurr<NTbuff
        buff(:, nsampcurr+1:NTbuff) = repmat(buff(:,nsampcurr), 1, NTbuff-nsampcurr);
    end
    if ops.GPU
        dataRAW = gpuArray(buff);
    else
        dataRAW = buff;
    end
    dataRAW = dataRAW';
    dataRAW = single(dataRAW);
    
    datr = filter(b1, a1, dataRAW);
    datr = flipud(datr);
    datr = filter(b1, a1, datr);
    datr = flipud(datr);
    
    bad=unique(bsxfun(@plus, find(sum(abs(datr)>sizeThresh,2)>chanThresh), -1:20));
    bad(bad<1)=[];
    bad(bad>size(datr,1))=[];
    
    datr(bad,:)=0;
    
    samples=gather_try(datr');
    fwrite(fidout, samples, 'int16');
    
    if ibatch<=Nbatch_buff
        isproc(ibatch) = 1;
    end
end

fclose(fid);
fclose(fidout);
ops.fbinary=fnameOut;
ops.filteredbinary=1;
