function [DecimalEvents, Timestamps] = OpenEphysEvents2Bpod(filename)
[data, pinChangeTimestamps, info] = load_open_ephys_data(filename);
Pos = find(info.eventId==1, 1, 'first');
BinaryEventCode = '0000000';
nPinChanges = length(pinChangeTimestamps)-Pos+1;
nTotalTimestamps = length(pinChangeTimestamps);
DecimalEvents = zeros(1,nPinChanges);
Timestamps = zeros(1,nPinChanges);
nEvents = 0;

while Pos <= nTotalTimestamps
    nPinsChanged = sum(pinChangeTimestamps == pinChangeTimestamps(Pos));
    for x = 1:nPinsChanged
         BinaryEventCode(8-(data(Pos))) = num2str(info.eventId(Pos));
%          BinaryEventCode(8-(data(Pos))) = num2str(info.eventId(Pos));

        Pos = Pos + 1;
    end
    nEvents = nEvents + 1;
    if (Pos <= nTotalTimestamps)
        DecimalEvents(nEvents) = bin2dec(BinaryEventCode);
        Timestamps(nEvents) = pinChangeTimestamps(Pos);
    else
        i = nPinChanges;
    end
end
DecimalEvents = DecimalEvents(1:nEvents-1);
Timestamps = Timestamps(1:nEvents-1);