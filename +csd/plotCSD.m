function plotCSD(stats, varargin)
% gets the event times of the current source density trials
% Inputs:
%   stats              [struct] - stats struct from csd.getCSD
% 
% ghs wrote it 2020

    ip = inputParser();
    ip.addParameter('overlayLFP', true)
    ip.addParameter('gamma', nan)
    ip.addParameter('numShanksPlot', nan)
    ip.parse(varargin{:});
    gamma = ip.Results.gamma;
    overlayLFP = ip.Results.overlayLFP;
    numShanksPlot = ip.Results.numShanksPlot;

    numShanks = size(stats.CSD, 3);
    

    if numShanks == 1
        if isfield(stats,'spectrogram')
            subplot(1,2,1)
        end
%         figure; clf;
%         imagesc(stats.time, stats.depth, stats.CSD-mean(stats.CSD(:))); axis ij
        imagesc(stats.time, stats.depth, stats.CSD); axis ij
        colormap(parula);
        hold on
        if overlayLFP
            plot(stats.time, bsxfun(@plus, stats.STA, stats.chDepths), 'Color', repmat(.1, 1, 3))
        end
        xlim(stats.time([1 end]))
%         plot(stats.latency, stats.sinkDepth(1), 'w*', 'Linewidth', 2)
%         plot(stats.latency, stats.sourceDepth(1), 'c*', 'Linewidth', 2)
        for e = 1:size(stats.reversalPointDepth{1}, 1)
            h = plot(stats.time([1 end]), [1; 1]*stats.reversalPointDepth{1}(e), 'r--', 'Linewidth', 2);
        end
        plot(stats.sinkTime{1}, stats.sinkDepth{1}, 'x')
        legend_labels = ['sink depth'];
        if ~any(isnan(stats.sinkDepth{1}))
%             yticks(stats.reversalPointDepth{shankInd})
            reversal_string = '';
            for e = 1:size(stats.reversalPointDepth{1}, 1)
                reversal_string = strcat(reversal_string, '\_', string(stats.reversalPointDepth{1}(e)));
            end
            legend_labels = [legend_labels strcat('reversal depth:  ', reversal_string)];

        end
        if isstruct(gamma)
            xLimits = get(gca,'XLim');
            scale = xLimits(2);
            plot(gamma.hgPower(:,1).*scale, stats.chDepths, 'Color', 'red', 'Linewidth', 1);
            legend_labels = [legend_labels 'high gamma'];
            plot(gamma.lgPower(:,1).*scale, stats.chDepths, 'Color', 'blue', 'Linewidth', 1);
            legend_labels = [legend_labels 'low gamma'];
        end
        l = legend(legend_labels);
        set(l,'FontSize',18);
        hold off
        colorbar
        
        if isfield(stats,'spectrogram')
%             figure(99); clf;
            subplot(1,2,2)
            imagesc(stats.time, stats.depth, stats.spectrogram(:,:,1)); hold on
            for e = 1:size(stats.reversalPointDepth{1}, 1)
                h = plot(stats.time([1 end]), [1; 1]*stats.reversalPointDepth{1}(e), 'r--', 'Linewidth', 2);
            end
            hold off
            colorbar
            xlabel('time')
            ylabel('depth')
        end
        
    elseif numShanks > 1
%         figure; clf;
        if ~isnan(numShanksPlot)
            numShanksLoop = min(numShanks, numShanksPlot);
        else
            numShanksLoop = numShanks;
        end
        for shankInd = 1:numShanksLoop
        %curShankInds = shankInd*lenShanks-lenShanks+1:shankInd*lenShanks;
        curCSD = stats.CSD(:,:,shankInd);
        if isfield(stats,'spectrogram')
            subplot(1,numShanks+2,shankInd)
        else
            subplot(1,numShanks,shankInd)
        end
%         imagesc(stats.time, stats.depth, curCSD-mean(curCSD(:))); axis ij
        imagesc(stats.time, stats.depth, curCSD); axis ij
        colormap(parula);
        hold on
        if overlayLFP
            plot(stats.time, bsxfun(@plus, stats.STA(:,:,shankInd), stats.chDepths), 'Color', repmat(.1, 1, 3))
        end
        xlim(stats.time([1 end]))
%         plot(stats.time([1 end]), stats.sinkDepth(shankInd)*[1 1], 'w--', 'Linewidth', 2)
%         plot(stats.time([1 end]), [1; 1]*stats.reversalPointDepth{shankInd}, 'r--', 'Linewidth', 2)
%         tmp = abs(stats.reversalPointDepth{1} - stats.sinkDepth(shankInd));
%         tmp = tmp + stats.sinkDepth(shankInd);
%         plot(stats.time([1 end]), [1; 1]*tmp, 'r--', 'Linewidth', 2)
        %                 axis ij
%         plot(stats.latency, stats.sinkDepth{shankInd}, 'w*', 'Linewidth', 2)
%         plot(stats.latency, stats.sourceDepth{shankInd}, 'c*', 'Linewidth', 2)
        plot(stats.sinkTime{shankInd}, stats.sinkDepth{shankInd}, 'x')
%         text(100,stats.sinkDepth{shankInd}+10,string(stats.sinkDepth{shankInd}))
        legend_labels = ['sink depth'];
        for e = 1:size(stats.reversalPointDepth{shankInd}, 1)
            h = plot(stats.time([1 end]), [1; 1]*stats.reversalPointDepth{shankInd}(e), 'r--', 'Linewidth', 2);
        end
        if ~any(isnan(stats.reversalPointDepth{shankInd}))
%             yticks(stats.reversalPointDepth{shankInd})
            reversal_string = '';
            for e = 1:size(stats.reversalPointDepth{shankInd}, 1)
                reversal_string = strcat(reversal_string, '\_', string(stats.reversalPointDepth{shankInd}(e)));
            end
            legend_labels = [legend_labels strcat('reversal depth:  ', reversal_string)];
            
        end
%         plot(stats.time([1 end]), [1; 1]*stats.reversalPointDepth{shankInd}, 'r--', 'Linewidth', 2)
        if isstruct(gamma)
            xLimits = get(gca,'XLim');
            scale = xLimits(2);
            plot(gamma.hgPower(:,shankInd).*scale, stats.chDepths, 'Color', 'red', 'Linewidth', 1);
            legend_labels = [legend_labels 'high gamma'];
            plot(gamma.lgPower(:,shankInd).*scale, stats.chDepths, 'Color', 'blue', 'Linewidth', 1);
            legend_labels = [legend_labels 'low gamma'];
        end
        l = legend(legend_labels);
        set(l,'FontSize',18);
        
        hold off
        colorbar
        end
        
        if isfield(stats,'spectrogram')
            for shankInd = 1:numShanks
%                 figure(99+shankInd); clf;
                subplot(1,numShanks+2,numShanks+shankInd)
                imagesc(stats.time, stats.depth, stats.spectrogram(:,:,shankInd)); hold on
                for e = 1:size(stats.reversalPointDepth{shankInd}, 1)
                    h = plot(stats.time([1 end]), [1; 1]*stats.reversalPointDepth{shankInd}(e), 'r--', 'Linewidth', 2);
                end
                hold off
                colorbar
                xlabel('time')
                ylabel('depth')
            end
        end
    end
end

