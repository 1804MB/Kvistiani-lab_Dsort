function [t, wv] = LoadTT_Intan(fn,records_to_get,record_units)
%LOADTT_INTAN   MClust loading engine for Intan data.
%   LOADTT_INTAN loads Intan data and fullfils all requirement for
%   Neuralynx loading engines (see MClust 3.5 documentation for details on
%   MClust loading engines).
%
%   Syntax:
%   [T, WV] = LOADTT_INTAN(FN,RECORDS_TO_GET,RECORD_UNITS)
%
%   Input arguments:
%       FN - file name with full path for Intan data to load
%       RECORDS_TO_GET, RECORD_UNITS - allows options for returning a
%           restricted data set; if not passed, the full data file is
%           returned;
%           RECORD_UNITS = 1 - return specific time stamps, given by
%               RECORDS_TO_GET
%           RECORD_UNITS = 2 - return data points at specific indices,
%               given by RECORDS_TO_GET
%           RECORD_UNITS = 3 - return data in a time stamp range, given by
%               RECORDS_TO_GET (2-elements vector)
%           RECORD_UNITS = 4 - return data in an index range, given by
%               RECORDS_TO_GET (2-elements vector)
%           RECORD_UNITS = 5 - return number of spikes as first output
%               argument; RECORDS_TO_GET should be empty
%
%   Output arguments:
%       T - time stamps in seconds; 1 x N, N = number of spikes (exeption:
%           RECORD_UNITS = 5, see above)
%       WV - waveforms; N X 4 X 30, N = number of spikes, 4 tetrode
%           channels, 30 time points for spikes
%
%   See also INTANDISC and LOADTT_NEURALYNXNT.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   9-May-2013

% Input argument check
error(nargchk(1,3,nargin))
switch nargin 
    case 1
        record_units = 0;
    case 2
        error('LoadTT_Intan:inputArg','Record_units argument should be provided.')
end

% Load data
TTdata = load(fn);
t = TTdata.TimeStamps;
wv = TTdata.WaveForms;

% Restrict data to the required range
switch record_units
    case 0   % return full data
        
    case 1   % return specific time stamps
        [jnk inxa] = intersect(t,records_to_get);  %#ok<*ASGLU> % indices for the time stamps
        t = t(inxa);
        wv = wv(inxa,:,:);
    case 2   % return specific data indices
        t = t(records_to_get);
        wv = wv(records_to_get,:,:);
    case 3   % select a time range
        if ~isequal(numel(records_to_get),2)
            error('LoadTT_Intan:inputArg','Input argument mismatch.')
        end
        [jnk inxa] = intersect(t,records_to_get);  % indices for the time stamps
        t = t(inxa(1):inxa(2));
        wv = wv(inxa(1):inxa(2),:,:);
    case 4   % select an index range
        if ~isequal(numel(records_to_get),2)
            error('LoadTT_Intan:inputArg','Input argument mismatch.')
        end
        t = t(records_to_get(1):records_to_get(2));
        wv = wv(records_to_get(1):records_to_get(2),:,:);   
    case 5   % special case: number of spikes
        if ~isempty(records_to_get)
            error('LoadTT_Intan:inputArg','Input argument mismatch.')
        end
        t = length(t);   % number of spikes
        wv = [];
    otherwise
        error('LoadTT_Intan:inputArg','Record_units should take an integer value from 1 to 5.')
end
if  record_units < 5
    t = t(:) * 1e4;   % enforce column vector form; use MClust time stamp convention
end