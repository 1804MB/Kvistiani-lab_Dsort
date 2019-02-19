% This function generates Trial Events file for Free Choice task using Bpod
% synchronization 
function MakeTrialEvents_FreeChoice(varargin)
% if nargin == 1
%     sessionpath = cellid2fnames(varargin{1},'sess');
% end
% cd(sessionpath)
[DecimalEvents, Timestamps] = OpenEphysEvents2Bpod('all_channels.events');
DecimalEvents(find(diff(Timestamps)<0.00005)+1) = [];
Timestamps(find(diff(Timestamps)<0.00005)+1) = [];
TrialInx = find(DecimalEvents == 0);

load('Jb.mat')
if length(TrialInx) ~= length(SessionData.TrialTypes)
    disp('DecimalEvents 0 do not match with SessionData TrialStart')
    TrialInx = TrialInx(1:length(SessionData.TrialTypes));
end
TrialStates = [];
inx = zeros(1,length(SessionData.RawData.OriginalStateNamesByNumber{1}));
val = cell(1,length(SessionData.RawData.OriginalStateNamesByNumber{1}));

TE = struct;
TE.CenterPortEntry = NaN(1,SessionData.nTrials);
TE.CenterPortExit = NaN(1,SessionData.nTrials);
TE.LeftPortEntry = NaN(1,SessionData.nTrials);
TE.LeftPortExit = NaN(1,SessionData.nTrials);
TE.RightPortEntry = NaN(1,SessionData.nTrials);
TE.RightPortExit = NaN(1,SessionData.nTrials);
TE.Reward = NaN(1,SessionData.nTrials);
TE.NoReward = NaN(1,SessionData.nTrials);
TE.TrialTypes = SessionData.TrialTypes;
TE.CenterPortLightOff = NaN(1,SessionData.nTrials);


for iTrial = 1:SessionData.nTrials
    if iTrial == 1
    TrialStates = [TrialStates SessionData.RawData.OriginalStateNamesByNumber{iTrial}(DecimalEvents(1:TrialInx(iTrial)-1)) 'EndTrial'];
    else
    TrialStates = [TrialStates SessionData.RawData.OriginalStateNamesByNumber{iTrial}(DecimalEvents(TrialInx(iTrial-1)+1:TrialInx(iTrial)-1)) 'EndTrial'];
    end
end
    
    TrialStatesByNumber = NaN(1,length(TrialStates));
    
    TrialStatesByNumber(strmatch('NewTrial',TrialStates)) = 1;
    TrialStatesByNumber(strmatch('InitialPoke',TrialStates)) = 2;
    TrialStatesByNumber(strmatch('WaitInCenter',TrialStates)) = 3;
    TrialStatesByNumber(strmatch('timeout',TrialStates)) = 4;
    TrialStatesByNumber(strmatch('WaitForPortExit',TrialStates)) = 5;
    TrialStatesByNumber(strmatch('WaitForReward',TrialStates)) = 6;
    TrialStatesByNumber(strmatch('FreeChoice',TrialStates)) = 7;
    TrialStatesByNumber(strmatch('NoRewardLeft',TrialStates)) = 8;
    TrialStatesByNumber(strmatch('NoRewardRight',TrialStates)) = 9;
    TrialStatesByNumber(strmatch('LeftReward',TrialStates)) = 10;
    TrialStatesByNumber(strmatch('Drink',TrialStates)) = 11;
    TrialStatesByNumber(strmatch('RightReward',TrialStates)) = 12;
    TrialStatesByNumber(strmatch('LeftInactivated',TrialStates)) = 13;
    TrialStatesByNumber(strmatch('RightInactivated',TrialStates)) = 14;
    TrialStatesByNumber(strmatch('DrinkingGrace',TrialStates)) = 15;
    TrialStatesByNumber(strmatch('EndTrial',TrialStates)) = 0;
    
    TrialEndInx = find(TrialStatesByNumber== 0);
    TrialStart = find(TrialStatesByNumber== 1);
	
    for iTrial = 1:length(TrialInx)   
	
    TrialStartInx(iTrial) = TrialStart(find(TrialStart < TrialEndInx(iTrial),1,'last'));  
    TrialInx = TrialStartInx(iTrial):TrialEndInx(iTrial);
	
    TE.CenterPortEntry(iTrial) = Timestamps(TrialInx(find(TrialStatesByNumber(TrialInx) == 1,1,'last'))) ;
    
	if ~isempty(TrialInx(TrialStatesByNumber(TrialInx) == 5))
        TE.CenterPortExit(iTrial) =  Timestamps(TrialInx(find(TrialStatesByNumber(TrialInx) == 5,1,'first'))) ;
        TE.CenterPortLightOff(iTrial) = Timestamps(TrialInx(find(TrialStatesByNumber(TrialInx) == 3,1,'first')));
    end
        
   
    if ~isempty(TrialInx(TrialStatesByNumber(TrialInx) == 8))
        TE.LeftPortEntry(iTrial) = Timestamps(TrialInx(TrialStatesByNumber(TrialInx) == 7)) ;
        TE.LeftPortExit(iTrial) = Timestamps(TrialInx(TrialStatesByNumber(TrialInx) == 8)) ;
        TE.NoReward(iTrial) = 1;
    end
    if ~isempty(TrialInx(TrialStatesByNumber(TrialInx) == 10))
        TE.LeftPortEntry(iTrial) = Timestamps(TrialInx(TrialStatesByNumber(TrialInx) == 7)) ;
        TE.LeftPortExit(iTrial) = Timestamps(TrialInx(find(TrialStatesByNumber(TrialInx) == 11,1,'first'))) ;
        TE.Reward(iTrial) = 1;
    end

    if ~isempty(TrialInx(TrialStatesByNumber(TrialInx) == 9))
        TE.RightPortEntry(iTrial) = Timestamps(TrialInx(TrialStatesByNumber(TrialInx) == 7)) ;
        TE.RightPortExit(iTrial) = Timestamps(TrialInx(TrialStatesByNumber(TrialInx) == 9)) ;        
        TE.NoReward(iTrial) = 1;
    end

    if ~isempty(TrialInx(TrialStatesByNumber(TrialInx) == 12))
        TE.RightPortEntry(iTrial) = Timestamps(TrialInx(TrialStatesByNumber(TrialInx) == 7)) ;
        TE.RightPortExit(iTrial) = Timestamps(TrialInx(find(TrialStatesByNumber(TrialInx) == 11,1,'first'))) ;     
        TE.Reward(iTrial) = 1;
    end
    

    clear TrialInx
    end 
TE.LeftPortOn = NaN(1,length(TE.TrialTypes));
TE.LeftPortOn(~isnan(TE.LeftPortEntry))=1;
TE.RightPortOn = NaN(1,length(TE.TrialTypes));
TE.RightPortOn(~isnan(TE.RightPortEntry))=1;

save('TrialEvents.mat','-struct','TE')
%  save TrialEvents TE