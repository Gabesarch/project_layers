%% WHY REMOVED SESSION
% 20, 21, 29, 43 no event times
% 23, 24, 25, 40, 41, 42 bad CSD
% 29, 31, 51 no CSD (nans)
% 37, 38 no import file
% 45 blob localization doesn't work

dataPath = getpref('FREEVIEWING', 'PROCESSED_DATA_DIR');

meta_file = fullfile(fileparts(which('addFreeViewingPaths')), 'Data', 'datasets.csv');
% meta_file = fullfile(dataPath, 'datasets.xls');

data = readtable(meta_file);

monkey = 'Logan';

sess_all = [1:19 22 26:28 30 32:36 39 44 46:50 53:57]; %30:35%22:28;%[1:19 21:30; %[40:42 44:47];
sess = [];
switch monkey
    case 'Logan'
        p = 1;
        for i = 1:length(sess_all)
            if strcmp(data{sess_all(i), 2}, monkey)
                sess(p) = sess_all(i);
                p = p+1;
            end
        end
    case 'Ellie'
        p = 1;
        for i = 1:length(sess_all)
            if strcmp(data{sess_all(i), 2}, monkey)
                sess(p) = sess_all(i);
                p = p+1;
            end
        end
    case 'both'
        sess = sess_all;
end

%% Get Average CSD

if strcmp(monkey, 'Ellie')
sess(sess==22) = [];
end

numTotPad = 50 %35 %230 %50 %35%230;
minDepth = -1500;

stats2 = [];
p = 1;

for i = sess
    
    disp(['Getting session ' num2str(i)])
    
    stats = csd.dataFactoryCSD(i, 'type', 'csd');
    
    gamma = csd.dataFactoryCSD(i, 'type', 'gamma');
    
    numShanks = size(stats.CSD, 3);
    
    for shankInd = 1:numShanks
        
        depDiff = abs(stats.depth(2)-stats.depth(1));
        
%         if isempty(stats.reversalPointDepth{shankInd})
%             stats.reversalPointDepth{shankInd} = stats.depth(1);
%         end
        
        baseInputDepth = gamma.lgInputLayerDepths(:,:,shankInd);
        baseInputDepth = baseInputDepth(1);
        
        %baseInputDepth = stats.reversalPointDepth{numShanks};
        
        depnew = stats.depth - baseInputDepth;
        
        amPad = round((minDepth - min(depnew))/(-depDiff));
        
        %depnew2 = [min(depnew)-depDiff*amPad:depDiff:min(depnew)-depDiff depnew max(depnew)+depDiff:depDiff:max(depnew)+depDiff*(numTotPad-amPad)];
        
        stats2(:,:,p) = [nan(numTotPad-amPad, size(stats.CSD,2)); stats.CSD(:,:,shankInd); nan(amPad, size(stats.CSD,2))];

        p = p+1;
    end
    
end

avgCSD = nanmean(stats2, 3);

depnew = linspace(-minDepth+(-depDiff)*size(avgCSD,1), -minDepth , size(avgCSD,1));



t = stats.time;

size(avgCSD)

figure; clf
imagesc(t, depnew, avgCSD); axis xy
colormap(jet)
colorbar
ylim([-700 1000])

%%

% figure(1);
depplot = (-24:-4)*50;
% p = polyfit(revPoint,troughs,1);
% f = polyval(p,revPoint);
c = lines;
c = c(1:4,:);
p = zeros(1,2);
dv = zeros(1,20);
for i = 1:length(conditionsT)
    [troughs, peaks, revPoint, T2Rev, P2Rev, TZPos, PZPos, ZAdj] = getBPPeaks(meta, conditionsP{i}, conditionsT{i}, 'CSDFlash.mat');
    
    figure(100);
    
    axes(ha(1))
    plot(revPoint,troughs,'ow', 'MarkerSize', 15, 'MarkerFaceColor', c(i,:)); hold on
    meanp1 = nanmean(T2Rev);
    stderror1= nanstd( T2Rev );
    
    inds2 = ~isnan(troughs);
    p = p + polyfit(revPoint(inds2),troughs(inds2),1);
    
    dv = dv+inds2;
    
    axes(ha(2))
    plot(revPoint,peaks,'ow', 'MarkerSize', 15, 'MarkerFaceColor', c(i,:)); hold on
%     scatter(revPoint,peaks, 'MarkerFaceColor', c(i,:)); hold on
    meanp1 = nanmean(P2Rev)
    stderror1= nanstd( P2Rev )
    
    
end

%% plot CSD with reversal and gamma points
% 25 26 27 41 52
% 56 57 do very well 
for i = sess; %[1:19 22:28 30 32:36 39:42 44:50 52:57] %44:50 %39:42  %32:36 %30 %22:28 %1:19 %[10:19 21:57]
    disp(['Getting session ' num2str(i)])
stats = csd.dataFactoryCSD(i, 'type', 'csd');
gamma = csd.dataFactoryCSD(i, 'type', 'gamma');
% figure(64);clf
% csd.plotCSD_BP(stats, gamma)
% pause
end

%% Plots of distance 

reversals = [];
troughs = [];
p = 1;
for i = sess
    
    disp(['Getting session ' num2str(i)])
    
    stats = csd.dataFactoryCSD(i, 'type', 'csd');
    
    gamma = csd.dataFactoryCSD(i, 'type', 'gamma');
    
    numShanks = size(stats.CSD, 3);
    
        for shankInd = 1:numShanks
            if isempty(stats.reversalPointDepth{shankInd})
                continue
            end
            reversals(p) = stats.reversalPointDepth{shankInd}(1);
            troughs(p) = gamma.lgTroughDepth(shankInd);
            
            p = p + 1;
        end
end

figure(10); clf

plot(reversals, troughs, 'o', 'Color','r'); hold on
refline([1 0])

% x = reversals;
% c = polyfit(x,troughs,1);
% y_est = polyval(c,x);
% plot(x,y_est,'r--','LineWidth',2)
xlabel('CSD Reversal Depth')
ylabel('Low Gamma Trough Depth')

%% Plots reversals/troughs that are wrong

reversals = [];
troughs = [];

reversalsBAD = [];
troughsBAD = [];
p = 1;
b = 1;
for i = sess
    
    disp(['Getting session ' num2str(i)])
    
    stats = csd.dataFactoryCSD(i, 'type', 'csd');
    
    gamma = csd.dataFactoryCSD(i, 'type', 'gamma');
    
    numShanks = size(stats.CSD, 3);
        
        
        for shankInd = 1:numShanks
            if isempty(stats.reversalPointDepth{shankInd})
                continue
            end
            reversals(p) = stats.reversalPointDepth{shankInd}(1);
            troughs(p) = gamma.lgTroughDepth(shankInd);
            
            if troughs(p)-reversals(p) > 300
                reversalsBAD(b) = reversals(p);
                troughsBAD(b) = troughs(p);
                b = b + 1;
                
                shankInd
                
                figure(74);clf
                csd.plotCSD_BP(stats, gamma) 
                pause
                
            end            
            p = p + 1;
        end
end

figure(10); clf

plot(reversals, troughs, 'o', 'Color','r'); hold on
refline([1 0])

% x = reversals;
% c = polyfit(x,troughs,1);
% y_est = polyval(c,x);
% plot(x,y_est,'r--','LineWidth',2)
xlabel('CSD Reversal Depth')
ylabel('Low Gamma Trough Depth')

    
    




