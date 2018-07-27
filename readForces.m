function [FxAll, FyAll, FzAll, TxAll, TyAll, TzAll, plateModAll, synchTrigAll, fs] = readForces(fullFilename,isTrig)

load(fullFilename);

fs = samplerate(1);
FxAll = data(datastart(1):dataend(1)); 
FyAll = data(datastart(2):dataend(2)); 
FzAll = data(datastart(3):dataend(3));
TxAll = data(datastart(4):dataend(4));
TyAll = data(datastart(5):dataend(5));
TzAll = data(datastart(6):dataend(6));

plateModAll = zeros(size(FxAll));
synchTrigAll = zeros(size(FxAll));
if isTrig
    plateModAll = data(datastart(7):dataend(7));
    synchTrigAll = data(datastart(8):dataend(8));
end
    
