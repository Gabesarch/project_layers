function stats = identifyMTSink(stats_sm1,stats_sm3)

    numShanks = size(stats_sm1.CSD, 3);
    ch0 = stats_sm1.depth;
    
    stats = stats_sm1;
    
    
%     if numShanks == 1
    for shankInd=1:numShanks
        CSD = stats_sm1.CSD(:,:,shankInd);
        
        time = stats_sm3.time;
        ix_time = time > 40 & time < 50;
        sm3_rvsl = stats_sm3.reversalPointDepth{shankInd}(1);
        depth = stats_sm3.depth;
        ix_depth = depth < sm3_rvsl & depth > sm3_rvsl - 400;
        [~,id] = min(reshape(CSD(ix_depth,ix_time), [], 1));
        
        % convert to indices
        [depthIndex,timeIndex] = ind2sub(size(CSD(:,ix_time)), id);
        
        time_ = time(ix_time);
        sinkTime = time_(timeIndex);
        
%         CSD = stats_sm1.CSD(:,:,1);
        % find reversal point
        CSD_ = CSD(:,ix_time);
        reversalPoints = findZeroCrossings(CSD_(:,timeIndex), -1);
%         indices_search = depth < depth(depthIndex) & depth > depth(find(ch0==sm3_rvsl)+1);
        indices_search = find(ix_depth);
        if length(indices_search)==0
            reversalPoints = [];
            disp('found no reversal')
        else
            reversalPoints = reversalPoints(reversalPoints>indices_search(1)-2 & reversalPoints<indices_search(end));
            reversalPoints = max(reversalPoints);
%             depthIndex = depthIndex+1;
        end
        
%         figure()
%         imagesc(CSD_(:,timeIndex))
        if isempty(reversalPoints)
            stats.reversalPointDepth{shankInd} = NaN;
            stats.sinkDepth{shankInd} = NaN;
            stats.sinkTime{shankInd} = NaN;
        else
            stats.reversalPointDepth{shankInd} = ch0(reversalPoints); %reversalPoints(1);
            stats.sinkDepth{shankInd} = ch0(depthIndex);
            stats.sinkTime{shankInd} = sinkTime;
        end
    end
    function i = findZeroCrossings(data, mode)
    %FINDZEROCROSSINGS Find zero crossing points.
    %   I = FINDZEROCROSSINGS(DATA,MODE) returns the indicies into the supplied
    %   DATA vector, corresponding to the zero crossings.
    %
    %   MODE specifies the type of crossing required:
    %     MODE < 0 - results in indicies for the -ve going zero crossings,
    %     MODE = 0 - results in indicies for ALL zero crossings (default), and
    %     MODE > 0 - results in indicies for the +ve going zero crossings.

    % $Id: findZeroCrossings.m,v 1.1 2008-07-21 23:31:50 shaunc Exp $

    if nargin < 2
        mode = 0;
    end

    [i,~,p] = find(data); % ignore zeros in the data vector

    switch sign(mode)
        case -1
            % find -ve going crossings
            ii = find(diff(sign(p))==-2);
        case 0
            % find all zero crossings
            ii = find(abs(diff(sign(p)))==2);
        case 1
            % find +ve going crossings
            ii = find(diff(sign(p))==2);
    end;

    i = round((i(ii)+i(ii+1))/2);
    end

    %         
    %     else
    %         
    %         
    %     end

    end

