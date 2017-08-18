function checkRawBinaryMatch(ops, NTbuff)
% checkRawBinaryMatch(ops, NTbuff)
if nargin < 2
    NTbuff = 1e3;
end

% read in binary data
fid = fopen(fullfile(ops.root, 'rawbinary.dat'));
buff = fread(fid, [ops.Nchan NTbuff], '*int16');

Dat = [];
for iCh = 1:ops.Nchan
    fl=dir(fullfile(ops.root, sprintf('100*CH%d.continuous', iCh)));
    fname = fullfile(ops.root, fl(1).name);
    [data, timestamps, info]=load_open_ephys_data_faster(fname);
    Dat = [Dat data(1:NTbuff)];
end

%%

figure(1); clf
subplot(1,2,1)
plot(Dat)
title('OE data')
xlabel('Sample #')
ylabel('mV')
Dat_=double(buff)'*ops.chHeaders(1).bitVolts;
subplot(1,2,2)
plot(Dat_)
MSE = mean( (Dat - Dat_).^2);
title('Data reconstructed from ops.rawbinary')
xlabel('Sample #')
ylabel('mV')
fprintf('Total MSE: %d\n', sum(MSE))
