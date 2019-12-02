function reqStats = RequestedInputs(officialStats,varargin)

list = varargin(1:end);
if (~isempty(list) && ~iscell(list{1}) && strcmpi(list{1},'all'))
    reqStats = 1:numel(officialStats);%officialStats;
elseif isempty(list)
    reqStats = 1:numel(officialStats);%officialStats;
else
    if (iscell(list{1}))
        list = list{1};
    end
    list = list(:);
    officialStatsL = lower(officialStats);
    reqStatsIdx = [];
    eid = sprintf('Images:%s:invalidMeasurement',mfilename);
    for k = 1:length(list)
      if (~ischar(list{k}))
        msg = sprintf('This measurement is not a string: "%d".', list{k});
        error(eid,'%s',msg);
      end
      idx = strmatch(lower(list{k}), officialStatsL);
      if (isempty(idx))
        msg = sprintf('Unknown measurement: "%s".', list{k});
        error(eid,'%s',msg);
      elseif (length(idx) > 1)
        msg = sprintf('Ambiguous measurement: "%s".', list{k});
        error(eid,'%s',msg);
      else
        reqStatsIdx = cat(1,reqStatsIdx,idx);
      end
    end
    reqStats = reqStatsIdx;%officialStats(reqStatsIdx);
end