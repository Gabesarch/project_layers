%% Import W

%% Add features to W struct

sess = [1:12 14:23 25 28:57]; %1:57;
dataPath = getpref('FREEVIEWING', 'PROCESSED_DATA_DIR');

meta_file = fullfile(fileparts(which('addFreeViewingPaths')), 'Data', 'datasets.csv');
% meta_file = fullfile(dataPath, 'datasets.xls');

data = readtable(meta_file);

error = [];
for sessionId = sess
    %try
    if isnumeric(sessionId)
        thisSession = data(sessionId,:);
    end
    
    S.processedFileName = [thisSession.Tag{1} '.mat'];
    
    fname_new = fullfile(dataPath, 'waveform', 'waveform_new', S.processedFileName);
    
    if ~exist(fname_new, 'file')
        disp('getting unprocessed Waveforms')
        fname = fullfile(dataPath, 'waveform', S.processedFileName);
    else
        disp('getting processed Waveforms')
        fname = fname_new;
    end
    
    
    fprintf('Loading [%s]\n', S.processedFileName)
    load(fname);
    fprintf('Done\n')
    
    
    
    %decide which waveforms to keep
    for i = 1:length(W)
        
        % normalize acutocorreleogram
        W(i).normACOR = (W(i).acor - min(W(i).acor)) / ( max(W(i).acor) - min(W(i).acor) );
        
        %get max time
        [~, mInd] = max(W(i).acor);
        W(i).tMaxAcor = W(i).tcor(mInd);
        
        %see if follows amplitude crit
        W(i).ampCrit = sqrt(sum(W(i).waveform(:,1).^2))/sqrt(sum(W(i).waveform(:,3).^2));
        
        % isi criterion
        if W(i).isiL < 0.73
            W(i).isi_keep = 1;
        else
            W(i).isi_keep = 0;
        end
        
        if W(i).ampCrit < 0.33
            W(i).amp_keep = 1;
        else
            W(i).amp_keep = 0;
        end
        
        
        % normalize by possion expectation
        % burst metric vs waveform duration (p2t)
        % find peak amplitude spike - find minimum amplitude spike - get
        % ratio between peak and minumum amplitude
        % cc = 1; sqrt(sum(W(cc).waveform(:,1).^2))/sqrt(sum(W(cc).waveform(:,3).^2))
        % Peak followed by trough
        % io.getwaveformstats - isiL
        % io.plotwaveform
        
    end
    
    W = rmfield(W,'peakInterp');
    
    save(fname_new, '-v7.3', 'W')
    %catch
    error = [error sessionId];
    %end
end

%% Add features to W struct part 2

sess = e3(4:end)

dataPath = getpref('FREEVIEWING', 'PROCESSED_DATA_DIR');

meta_file = fullfile(fileparts(which('addFreeViewingPaths')), 'Data', 'datasets.csv');
% meta_file = fullfile(dataPath, 'datasets.xls');

data = readtable(meta_file);

error = [];
for sessionId = sess
    %try
        disp(['processing session: ' num2str(sessionId)])
        if isnumeric(sessionId)
            thisSession = data(sessionId,:);
        end
        
        S.processedFileName = [thisSession.Tag{1} '.mat'];
        
        disp('getting Waveforms')
        
        fname_new = fullfile(dataPath, 'waveform', 'waveform_new', S.processedFileName);
        
        fprintf('Loading [%s]\n', S.processedFileName)
        
        load(fname_new);
        fprintf('Done\n')
        
        Exp = io.dataFactoryGratingSubspace(sessionId);
        
        % trial starts and stops
        tstart = Exp.ptb2Ephys(cellfun(@(x) x.STARTCLOCKTIME, Exp.D(:)));
        tstop = Exp.ptb2Ephys(cellfun(@(x) x.ENDCLOCKTIME, Exp.D(:)));
        
        % inter-trial intervals
        epochs = [tstop(1:end-1) tstart(2:end)];
         
        disp('Getting waveform stats')
        % get waveform statistics over valid epochs
        W_ = io.get_waveform_stats(Exp.osp, 'validEpochs', epochs);
        cids = [W.cid];
        cidsW_ = [W_.cid];
        for i = 1:length(cids)
            ind = find(cidsW_==cids(i));
            W(i).isiRate = W_(ind).isiRate;
            W(i).localityIdx = W_(ind).localityIdx;
            W(i).BRIJake = W_(ind).BRI;
            W(i).normrate = W_(ind).normrate;
        end
        
        save(fname_new, '-v7.3', 'W')
    %catch
        %error = [error sessionId];
    %end
end

%% Get W_all

sess = [1:12 14:23 25 28:36 39:57]; %[1:12 14:23 25 28:57];
dataPath = getpref('FREEVIEWING', 'PROCESSED_DATA_DIR');

meta_file = fullfile(fileparts(which('addFreeViewingPaths')), 'Data', 'datasets.csv');
% meta_file = fullfile(dataPath, 'datasets.xls');

data = readtable(meta_file);

W_all = [];
for sessionId = sess
    
    if isnumeric(sessionId)
        thisSession = data(sessionId,:);
    end
    
    S.processedFileName = [thisSession.Tag{1} '.mat'];
    
    fname_new = fullfile(dataPath, 'waveform', 'waveform_new', S.processedFileName);
    
    load(fname_new);
    
    W_all = [W_all;W];
    
end



%% Get criterion values
figure(1); clf
subplot(1,2,1)
isiL = [W_all.isiL];
isiL(isiL>1) = nan;
histogram(isiL, 100); hold on
plot(0.73*[1 1], ylim, 'r')
title('Line noise in ISI')
legend('hist', 'line noise threshold')


isiL_keep = isiL(isiL<0.73);
% 0.73


subplot(1,2,2)
ampCrit = [W_all.ampCrit];
ampCrit = ampCrit(ampCrit<1);
histogram(ampCrit,50); hold on
plot(0.33*[1 1], ylim, 'r')
title('Waveform localization metric')
legend('hist', 'Waveform localization threshold')
% 0.33

%% get layer info for each cell

% stopped at 17
sess = 17:57

dataPath = getpref('FREEVIEWING', 'PROCESSED_DATA_DIR');

meta_file = fullfile(fileparts(which('addFreeViewingPaths')), 'Data', 'datasets.csv');
% meta_file = fullfile(dataPath, 'datasets.xls');

data = readtable(meta_file);

error = [];
for sessionId = sess
    try
        if isnumeric(sessionId)
            thisSession = data(sessionId,:);
        end
        
        S.processedFileName = [thisSession.Tag{1} '.mat'];
        
        disp('getting Waveforms')
        
        fname_new = fullfile(dataPath, 'waveform', 'waveform_new', S.processedFileName);
        
        fprintf('Loading [%s]\n', S.processedFileName)
        
        load(fname_new);
        fprintf('Done\n')
        
        disp('getting gamma')
        fname = fullfile(dataPath, 'gamma', S.processedFileName);
        
        fprintf('Loading [%s]\n', S.processedFileName)
        
        gamma = load(fname);
        fprintf('Done\n')
        
        fname_new = fullfile(dataPath, 'waveform', 'waveform_new', S.processedFileName);
        
        load(fname_new);
        
        disp(['Getting spikes: session ' num2str(sessionId)])
        
        spfilename = [thisSession.Tag{1} '_kilo.mat'];
        spfname = fullfile(dataPath, 'spikes', spfilename);
        sp = load(spfname);
        sp = sp.sp; %CHANGE
        
        for i = 1:length(W)
            W(i).spikes = sp.st(sp.clu==W(i).cid);
            
            [W(i).tcor4,W(i).acor4,W(i).pcor4,W(i).ncor4, W(i).whichonesempty4] = comp_autocor_fast(W(i).spikes,0.0004,0.040);
            [W(i).tcor_smooth4,W(i).acor_smooth4,W(i).pcor_smooth4,W(i).ncor_smooth4, ~,W(i).BI,W(i).PK,W(i).LENO] = comp_autocor_fast_JM(W(i).spikes,0.0004,0.040);
            %[W(i).tcor2,W(i).acor2,W(i).pcor2,W(i).ncor2, ~] = comp_autocor_fast(W(i).spikes,0.0002,0.040);
            %[W(i).tcor_smooth2,W(i).acor_smooth2,W(i).pcor_smooth2,W(i).ncor_smooth2, ~,W(i).BI2,W(i).PK2,W(i).LENO2] = comp_autocor_fast_JM(W(i).spikes,0.0002,0.040);
            
            % Decide which layer unit is in
            
            
            W(i).inputDepths = gamma.lgInputLayerDepths(:,:,1);
            if W(i).depth > max(W(i).inputDepths)
                W(i).layer = 3; % cell is in deep layer
            elseif W(i).depth < max(W(i).inputDepths) && W(i).depth > min(W(i).inputDepths)
                W(i).layer = 2; % cell is in input layer
            elseif W(i).depth < min(W(i).inputDepths)
                W(i).layer = 1; % cell is in superficial layer
            end
            
            
        end
        
        save(fname_new, '-v7.3', 'W')
    catch
        error = [error sessionId];
    end
end

%% Add in autocorr and spikes to W struct

sess = 1:57

dataPath = getpref('FREEVIEWING', 'PROCESSED_DATA_DIR');

meta_file = fullfile(fileparts(which('addFreeViewingPaths')), 'Data', 'datasets.csv');
% meta_file = fullfile(dataPath, 'datasets.xls');

data = readtable(meta_file);

error = [];
for sessionId = sess
    try
        if isnumeric(sessionId)
            thisSession = data(sessionId,:);
        end
        
        S.processedFileName = [thisSession.Tag{1} '.mat'];
        
        disp('getting Waveforms')
        
        fname_new = fullfile(dataPath, 'waveform', 'waveform_new', S.processedFileName);
        
        fprintf('Loading [%s]\n', S.processedFileName)
        
        load(fname_new);
        fprintf('Done\n')
        
        disp(['Getting spikes: session ' num2str(sessionId)])
        
        spfilename = [thisSession.Tag{1} '_kilo.mat'];
        spfname = fullfile(dataPath, 'spikes', spfilename);
        sp = load(spfname);
        sp = sp.sp; %CHANGE
        
        %         Keep1 = [W.isi_keep];
        %         Keep2 = [W.amp_keep];
        %         Keep = Keep1 + Keep2;
        %
        %         tmpW = W(find(Keep == 2));
        
        for i = 1:length(W)
            W(i).spikes = sp.st(sp.clu==W(i).cid);
        end
        
        [W(i).tcor4,W(i).acor4,W(i).pcor4,W(i).ncor4, W(i).whichonesempty4] = comp_autocor_fast(W(i).spikes,0.0004,0.040);
        [W(i).tcor_smooth4,W(i).acor_smooth4,W(i).pcor_smooth4,W(i).ncor_smooth4, ~,W(i).BI,W(i).PK,W(i).LENO] = comp_autocor_fast(W(i).spikes,0.0004,0.040);
        [W(i).tcor4,W(i).acor4,W(i).pcor4,W(i).ncor4, W(i).whichonesempty4] = comp_autocor_fast(W(i).spikes,0.0004,0.040);
        [W(i).tcor_smooth4,W(i).acor_smooth4,W(i).pcor_smooth4,W(i).ncor_smooth4, ~,W(i).BI,W(i).PK,W(i).LENO] = comp_autocor_fast(W(i).spikes,0.0004,0.040);
        save(fname_new, 'W')
    catch
        error = [error sessionId];
    end
end

%% Add in autocorr and spikes to W struct

sess = 1:57

dataPath = getpref('FREEVIEWING', 'PROCESSED_DATA_DIR');

meta_file = fullfile(fileparts(which('addFreeViewingPaths')), 'Data', 'datasets.csv');
% meta_file = fullfile(dataPath, 'datasets.xls');

data = readtable(meta_file);

error = [];
for sessionId = sess
    try
        if isnumeric(sessionId)
            thisSession = data(sessionId,:);
        end 
        
        disp(['Processing session: ' num2str(sessionId)])
        
        S.processedFileName = [thisSession.Tag{1} '.mat'];
        
        disp('getting Waveforms')
        
        fname_new = fullfile(dataPath, 'waveform', 'waveform_new', S.processedFileName);
        
        fprintf('Loading [%s]\n', S.processedFileName)
        
        load(fname_new);
        fprintf('Done\n')
        
        for i = 1:length(W)
            [W(i).tcor_logsmooth,W(i).acor_logsmooth,W(i).pcor_logsmooth,W(i).ncor_logsmooth,W(i).BI_logsmooth,W(i).PK_logsmooth,W(i).LENO_logsmooth] = comp_autocor_fit(W(i).spikes,0.0004,0.040);
        end
        
        save(fname_new, 'W')
    catch
        error = [error sessionId];
    end
end

%% Add in spikes to W struct

sess = 1:57

dataPath = getpref('FREEVIEWING', 'PROCESSED_DATA_DIR');

meta_file = fullfile(fileparts(which('addFreeViewingPaths')), 'Data', 'datasets.csv');
% meta_file = fullfile(dataPath, 'datasets.xls');

data = readtable(meta_file);

error = [];
for sessionId = sess
    try
        if isnumeric(sessionId)
            thisSession = data(sessionId,:);
        end
        
        S.processedFileName = [thisSession.Tag{1} '.mat'];
        
        disp('getting Waveforms')
        
        fname_new = fullfile(dataPath, 'waveform', 'waveform_new', S.processedFileName);
        
        fprintf('Loading [%s]\n', S.processedFileName)
        
        load(fname_new);
        fprintf('Done\n')
        
        Exp = io.dataFactoryGratingSubspace(i);
        
        sp = Exp.osp;
        
        sp.st(sp.clu==p)
    catch
        error = [error sessionId];
    end
end



%% Get W_keep - START HERE IF ALREADY GENERATED FILES

sess = [1:12 14:23 25 28:36 39:57];% 14:23 25 28:57];
dataPath = getpref('FREEVIEWING', 'PROCESSED_DATA_DIR');

meta_file = fullfile(fileparts(which('addFreeViewingPaths')), 'Data', 'datasets.csv');
% meta_file = fullfile(dataPath, 'datasets.xls');

data = readtable(meta_file);

W_keep = [];
for sessionId = sess
    
    if isnumeric(sessionId)
        thisSession = data(sessionId,:);
    end
        
    S.processedFileName = [thisSession.Tag{1} '.mat'];
    
    fname_new = fullfile(dataPath, 'waveform', 'waveform_new', S.processedFileName);
    
    load(fname_new);
    
     disp(['Getting waveforms: session ' num2str(sessionId)])
%     
%     spfilename = [thisSession.Tag{1} '_kilo.mat'];
%     spfname = fullfile(dataPath, 'spikes', spfilename);
%     sp = load(spfname);
%     sp = sp.sp; %CHANGE
    isiR = arrayfun(@(x) x.isiRate, W); % refractory rate / expected rate
    localIdx = arrayfun(@(x) x.localityIdx, W); % wf amplitude ch-2 / ch center
    ix = localIdx < .5 & isiR<1;% & p2t > 0;
    
%     Keep1 = [W.isi_keep];
%     Keep2 = [W.amp_keep];
%     Keep = Keep1 + Keep2;
    
    tmpW = W(ix);
    
%     for i = 1:length(tmpW)
%         tmpW(i).spikes = sp.st(sp.clu==tmpW(i).cid);
%     end
    
    W_keep = [W_keep; tmpW];
    
end

%% Get more accurate peak to trough time & only keep with positive p2t

plotIt = false;

for i = 1:length(W_keep)
    xq = linspace(W_keep(i).wavelags(1), W_keep(i).wavelags(end), 500);
    vq = interp1(W_keep(i).wavelags,W_keep(i).waveform(:,3),xq, 'spline');
    vqabs = abs(vq);
    
    % [~,ind] = max(vq);
    % W(i).peakInterp = xq(ind);
    
    [~,P] = islocalmax(vqabs);
    
    %I = find(TF);
    
    [max1, ind1] = max(P);
    %P(1:ind1) = -1;
    %W_keep(i).peak1 = xq(ind1);
    
    % second peak must be opposite sign of first peak
    if sign(vq(ind1))== -1
        vqsign = vq;
        vqsign(vqsign<0) = 0;
    elseif sign(vq(ind1))== 1
        vqsign = vq;
        vqsign(vqsign>0) = 0;
    end
    
    [~,P] = islocalmax(abs(vqsign));
    
    P(1:ind1) = 0; % must find second peak after first peak
    
    [max2, ind2] = max(P);
    if ind1>ind2
        W_keep(i).peak2time = xq(ind1);
        W_keep(i).peak1time = xq(ind2);
        W_keep(i).peak2val = vq(ind1);
        W_keep(i).peak1val = vq(ind2);
    elseif ind2>ind1
        W_keep(i).peak2time = xq(ind2);
        W_keep(i).peak1time = xq(ind1);
        W_keep(i).peak2val = vq(ind2);
        W_keep(i).peak1val = vq(ind1);
    end
    
    W_keep(i).t2p_interp = [W_keep(i).peak2time-W_keep(i).peak1time];
    W_keep(i).peakvaldiff = [W_keep(i).peak2val-W_keep(i).peak1val];
    
    if plotIt
        figure(1); clf
        plot(W_keep(i).wavelags, W_keep(i).waveform(:,3)); hold on
        %plot(xq, vq);
        plot(W_keep(i).peak1*[1 1], ylim, 'r')
        plot(W_keep(i).peak2*[1 1], ylim, 'r')
        pause
    end
end

%t2p_all = [W_keep.t2p_interp];
t2p_all = arrayfun(@(x) x.t2p_interp, W_keep); % refractory rate / expected rate

peakvaldiff = arrayfun(@(x) x.peakvaldiff, W_keep); % refractory rate / expected rate

ix = peakvaldiff > 0;

W_keep = W_keep(ix);



% figure(30); clf
% plot(xq, vq ); hold on;
% plot(xq(ind)*[1 1], ylim, 'r');
% plot(W(i).wavelags, W(i).waveform(:,3)); hold off
% pause

%% Get peak2trough distribution - ALL

%peaksT = [W_all.peaktime];
%troughsT = [W_all.troughtime];
%thresh = 0.0004063;
%thresh = 0.000503;
thresh = 0.0003552;


t2p_all = [W_keep.t2p_interp];
figure(29); clf
histogram(t2p_all*1000, 72); hold on
plot(thresh*1000*[1 1], ylim, 'r')
title('Peak-to-Trough Distribution')
legend('t2p', 'NW-BW threshold')
xlabel('trough-to-peak time (ms)')

% 0.0004063

%% Get NW cells and BW cells

t2p_all = [W_keep.t2p_interp];

%thresh = 0.000503;

thresh = thresh;

corr_peak_times = [W_keep.tMaxAcor];

%NW
NWcells = find(t2p_all<thresh);
[~, I] = sort(corr_peak_times);
I = intersect( I, NWcells);
[~, I2] = sort(corr_peak_times(I));
I_NW = I(I2);

% sort indices by peak NW time
W_NW = [];
for i = I_NW
    W_NW = [W_NW; W_keep(i)];
end

%BW
BWcells = find(t2p_all>=thresh);
[~, I] = sort(corr_peak_times);
I = intersect( I, BWcells);
[~, I2] = sort(corr_peak_times(I));
I_BW = I(I2);

% sort indices by peak NW time
W_BW = [];
for i = I_BW
    W_BW = [W_BW; W_keep(i)];
end

%% Plot each autocorrelogram one at a time
for i = 1:length(W_NW)
    figure(1); clf;
    plot(W_NW(i).tcor,W_NW(i).acor,'k.-')
    pause
end

%% Plot NW and BW autocorr
tcor = W_keep(1).tcor;

figure(3); clf
subplot(1,2,1)
acor_NW = [];
for i = 1:length(W_NW)
    acor_NW = [acor_NW; W_NW(i).normACOR];
end
imagesc(tcor, 1:length(acor_NW), acor_NW)
title('NW cells only')

subplot(1,2,2)
acor_BW = [];
for i = 1:length(W_BW)
    acor_BW = [acor_BW; W_BW(i).normACOR];
end
imagesc(tcor, 1:length(acor_BW), acor_BW)
title('BW cells only')

%%

%% Compute burstiness/refractoriness index
BRI_NW = [];
start = 1;
last = 4;
for i = 1:length(W_NW)
    BRI_NW(i) = mean(W_NW(i).ncor((W_NW(i).tcor>=start)&(W_NW(i).tcor<=last)));
    
end
I = find(BRI_NW<20);

BRI_NW = BRI_NW(I);

figure(31); clf
histogram(BRI_NW, 200)

p2t_NW = [W_NW.t2p_interp];
p2t_NW = p2t_NW(I);

W_tmp = W_NW(I);
figure(34); clf
plot(p2t_NW, BRI_NW, '.', 'MarkerSize', 8); hold on
%set(gca, 'YScale', 'log')
xlabel('Waveform duration (s)')
ylabel('Peak of Normalized Autocorrelation')

%% Compute burstiness/refractoriness index - BY LAYER
BRI_NW = [];
start = 1;
last = 4;
for i = 1:length(W_NW)
    BRI_NW(i) = mean(W_NW(i).ncor((W_NW(i).tcor>=start)&(W_NW(i).tcor<=last)));
    
end
I = find(BRI_NW<20);

BRI_NW = BRI_NW(I);

p2t_NW = [W_NW.t2p_interp];
p2t_NW = p2t_NW(I);

W_tmp = W_NW(I);
figure(35); clf

for i = 1:3
    plot(p2t_NW([W_tmp.layer]==i), BRI_NW([W_tmp.layer]==i), '.', 'MarkerSize', 8); hold on
end
%set(gca, 'YScale', 'log')
legend('Superficial Layer', 'Input Layer', 'Deep Layer')
xlabel('Waveform duration (s)')
ylabel('Peak of Normalized Autocorrelation')

% BRI_NW = [];
% for i = 1:length(W_NW)
%     BRI_NW(i) = W_NW(i).tMaxAcor/W_NW(i).pcor(1);
% end
% %BRI_NW(BRI_NW>10) = nan;
% histogram(BRI_NW, 200)




%% AutoCorr - sorted by max autocorr peak time

for i = 1:length(W_keep)
    W_keep(i).normACOR = (W_keep(i).ncor - min(W_keep(i).ncor)) / ( max(W_keep(i).ncor) - min(W_keep(i).ncor) );
end
corr_peak_times = [W_keep.tMaxAcor];
[~, I] = sort(corr_peak_times);

acor = [];
for i = I
    acor = [acor; W_keep(i).normACOR];
end
tcor = W_keep(1).tcor;

imagesc(tcor, 1:length(W_keep), acor)


%% Waveform sanity check
waveformsNW = [];
for i = 1:length(W_NW)
    wav = W_NW(i).waveform(:,3);
    [~, ind] = max(abs(wav));
    shiftnum = 31-ind;
    wav = circshift(wav,shiftnum);
    waveformsNW = [waveformsNW wav];
end

waveformsBW = [];
for i = 1:length(W_BW)
    wav = W_BW(i).waveform(:,3);
    [~, ind] = max(abs(wav));
    shiftnum = 31-ind;
    wav = circshift(wav,shiftnum);
    waveformsBW = [waveformsBW wav];
end

figure(45); clf
plot(waveformsNW, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 0.001); hold on
plot(waveformsBW, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 0.001);


%%
W4 = [];
emptyones = [];
sess = [1:12 14:23 25 28:57];
for i = 1:length(W_keep)
    
    %[tcor1,acor1,pcor1,ncor1] = comp_autocor_fast(W_keep(i).spikes,0.0001,0.040);
    %[tcor2,acor2,pcor2,ncor2] = comp_autocor_fast(W_keep(i).spikes,0.0002,0.040);
    %[tcor3,acor3,pcor3,ncor3] = comp_autocor_fast(W_keep(i).spikes,0.0003,0.040);
    %         disp(size(W_keep(i).spikes))
    %         W_keep(i).size = size(W_keep(i).spikes);
    %         [W_keep(i).tcor4,W_keep(i).acor4,W_keep(i).pcor4,W_keep(i).ncor4, W_keep(i).whichonesempty] = comp_autocor_fast(W_keep(i).spikes,0.0004,0.040);
    [W_keep(i).tcor4,W_keep(i).acor4,W_keep(i).pcor4,W_keep(i).ncor4] = comp_autocor_fast(W_keep(i).spikes,0.0004,0.040);
    
    
    %         if W_keep(i).whichonesempty
    %             emptyones = [emptyones i];
    %         end
    
end

save('W_keep', 'W_keep')


%%


%% Compute burstiness/refractoriness index
BRI_NW = [];
start = 1;
last = 4;
for i = 1:length(W_NW)
    BRI_NW(i) = mean(W_NW(i).ncor_smooth4((W_NW(i).tcor_smooth4>=start)&(W_NW(i).tcor_smooth4<=last)));
    %     BRI_NW(i) = max(W_NW(i).ncor_smooth4);%((W_NW(i).tcor_smooth4>=start)&(W_NW(i).tcor_smooth4<=last)));
    %     plot(W_NW(i).ncor_smooth4)
    %     pause
end
I = find(BRI_NW<20);

BRI_NW = BRI_NW(I);

figure(31); clf
histogram(BRI_NW, 200)

p2t_NW = [W_NW.t2p_interp];
p2t_NW = p2t_NW(I);

BI_NW = [W_NW.BI];
BI_NW = BI_NW(I);

W_tmp = W_NW(I);
figure(34); clf
% plot(p2t_NW, BRI_NW, '.', 'MarkerSize', 8); hold on
% plot(BI_NW, BRI_NW, '.', 'MarkerSize', 8); hold on
plot(BI_NW, p2t_NW, '.', 'MarkerSize', 8); hold on


%set(gca, 'YScale', 'log')
xlabel('Waveform duration (s)')
ylabel('Peak of Normalized Autocorrelation')

%% Compute burstiness/refractoriness index - BY LAYER
BRI_NW = [];
start = 1;
last = 4;
for i = 1:length(W_NW)
    BRI_NW(i) = mean(W_NW(i).ncor4((W_NW(i).tcor4>=start)&(W_NW(i).tcor4<=last)));
    %     BRI_NW(i) = max(W_NW(i).ncor4((W_NW(i).tcor4>=start)&(W_NW(i).tcor4<=last)));
    
    
end
I = find(BRI_NW<20);

BRI_NW = BRI_NW(I);

p2t_NW = [W_NW.t2p_interp];
p2t_NW = p2t_NW(I);

W_tmp = W_NW(I);
figure(35); clf

for i = 1:3
    plot(p2t_NW([W_tmp.layer]==i), BRI_NW([W_tmp.layer]==i), '.', 'MarkerSize', 8); hold on
end
%set(gca, 'YScale', 'log')
legend('Superficial Layer', 'Input Layer', 'Deep Layer')
xlabel('Waveform duration (s)')
ylabel('Peak of Normalized Autocorrelation')

%% Compute burstiness/refractoriness index - BY LAYER - graded colors
BRI_NW = [];
start = 1;
last = 4;
meth = 'max'; % method for computing BRI

%var1 = 'BRI'
%var2 = 'p2t'

var2 = 'BRI';
var1 = 'p2t';

var1tit = var1;
var2tit = var2;
tit = ['var1=' var1 ',var2=' var2];


for i = 1:length(W_NW)
    switch meth
        case '1to4'
            BRI_NW(i) = mean(W_NW(i).ncor_smooth4((W_NW(i).tcor_smooth4>=start)&(W_NW(i).tcor_smooth4<=last)));
        case 'max'
            [~,BRI_NW(i)] = max(W_NW(i).ncor_smooth4);
        case '1-3/4-6'
            BRI_NW(i) = mean(W_NW(i).ncor_smooth4((W_NW(i).tcor_smooth4>=1)&(W_NW(i).tcor_smooth4<=3)))/mean(W_NW(i).ncor_smooth4((W_NW(i).tcor_smooth4>=4)&(W_NW(i).tcor_smooth4<=6)));
        case 'acor1to4'
            BRI_NW(i) = mean(W_NW(i).acor_smooth4((W_NW(i).tcor_smooth4>=start)&(W_NW(i).tcor_smooth4<=last)));
        case 'acormax'
            BRI_NW(i) = max(W_NW(i).acor_smooth4);
        case 'acor1-3/4-6'
            BRI_NW(i) = mean(W_NW(i).acor_smooth4((W_NW(i).tcor_smooth4>=1)&(W_NW(i).tcor_smooth4<=3)))/mean(W_NW(i).acor_smooth4((W_NW(i).tcor_smooth4>=4)&(W_NW(i).tcor_smooth4<=6)));
    end
    %     BRI_NW(i) = max(W_NW(i).ncor4((W_NW(i).tcor4>=start)&(W_NW(i).tcor4<=last)));
    
    
end

dist2input = [W_NW.inputDepths]; %input layer depths
dist2input = dist2input(1:2:end);
dist2input = dist2input - [W_NW.depth];% get distance to reversal - negative number is below reversal

switch var1
    case 'BRI'
        var1 = BRI_NW;
        I = find(BRI_NW<40);
        tit = [tit ',method=' meth];
    case 'p2t'
        var1 = [W_NW.t2p_interp];
        I = 1:length(var1);
    case 'BI'
        var1 = [W_NW.BI];
        I = 1:length(var1);
    case 'laminar'
        var1 = dist2input;
        I = 1:length(var1);
end

switch var2
    case 'BRI'
        var2 = BRI_NW;
        I = find(BRI_NW<40);
        tit = [tit ',method=' meth];
    case 'p2t'
        var2 = [W_NW.t2p_interp];
    case 'BI'
        var2 = [W_NW.BI];
    case 'laminar'
        dist2input = [W_NW.inputDepths]; %input layer depths
        dist2input = dist2input(1:2:end);
        dist2input = dist2input - [W_NW.depth];% get distance to reversal - negative number is below reversal
        var2 = dist2input;
end

length(I)
var1 = var1(I);
var2 = var2(I);
dist2input = dist2input(I);


% I = find(BRI_NW<20);
%
% BRI_NW = BRI_NW(I);
%
% p2t_NW = [W_NW.t2p_interp];
% p2t_NW = p2t_NW(I);
%
% W_tmp = W_NW(I);
% figure(35); clf
%
% dist2input = [W_NW.inputDepths]; %input layer depths
% dist2input = dist2input(1:2:end);
% dist2input = dist2input - [W_NW.depth];% get distance to reversal - negative number is below reversal
% dist2input = dist2input(I);

figure(87); clf
scatter(var2, var1, 36, dist2input, 'filled'); hold on
set(gca,'yscale','log')
set(gca,'xscale','log')
colormap jet
colorbar

%set(gca, 'YScale', 'log')
%legend('Superficial Layer', 'Input Layer', 'Deep Layer')
title(tit)
xlabel(var2tit)
ylabel(var1tit)
hold off

%%

%% Compute burstiness/refractoriness index - BY LAYER - graded colors
BRI_NW = [];
start = 1;
last = 4;
meth = '1to4'; % method for computing BRI

%var1 = 'BRI'
%var2 = 'p2t'

var2 = 'BI';
var1 = 'BRI';

var1tit = var1;
var2tit = var2;
tit = ['var1=' var1 ',var2=' var2];


for i = 1:length(W_NW)
    switch meth
        case '1to4'
            BRI_NW(i) = mean(W_NW(i).ncor_logsmooth((W_NW(i).tcor_logsmooth>=start)&(W_NW(i).tcor_smooth4<=last)));
        case 'max'
            [~,BRI_NW(i)] = max(W_NW(i).ncor_logsmooth);
        case '1-3/4-6'
            BRI_NW(i) = mean(W_NW(i).ncor_logsmooth((W_NW(i).tcor_logsmooth>=1)&(W_NW(i).tcor_logsmooth<=3)))/mean(W_NW(i).ncor_logsmooth((W_NW(i).tcor_logsmooth>=4)&(W_NW(i).tcor_logsmooth<=6)));
        case 'acor1to4'
            BRI_NW(i) = mean(W_NW(i).acor_logsmooth((W_NW(i).tcor_loglogsmooth>=start)&(W_NW(i).tcor_logsmooth<=last)));
        case 'acormax'
            BRI_NW(i) = max(W_NW(i).acor_logsmooth);
        case 'acor1-3/4-6'
            BRI_NW(i) = mean(W_NW(i).acor_logsmooth((W_NW(i).tcor_loglogsmooth>=1)&(W_NW(i).tcor_logsmooth<=3)))/mean(W_NW(i).acor_logsmooth((W_NW(i).tcor_logsmooth>=4)&(W_NW(i).tcor_logsmooth<=6)));
    end
    %     BRI_NW(i) = max(W_NW(i).ncor4((W_NW(i).tcor4>=start)&(W_NW(i).tcor4<=last)));
    
    
end

dist2input = [W_NW.inputDepths]; %input layer depths
dist2input = dist2input(1:2:end);
dist2input = dist2input - [W_NW.depth];% get distance to reversal - negative number is below reversal

switch var1
    case 'BRI'
        var1 = BRI_NW;
        I = find(BRI_NW<40);
        tit = [tit ',method=' meth];
    case 'p2t'
        var1 = [W_NW.t2p_interp];
        I = 1:length(var1);
    case 'BI'
        var1 = [W_NW.BI_logsmooth];
        I = 1:length(var1);
    case 'laminar'
        var1 = dist2input;
        I = 1:length(var1);
end

switch var2
    case 'BRI'
        var2 = BRI_NW;
        I = find(BRI_NW<40);
        tit = [tit ',method=' meth];
    case 'p2t'
        var2 = [W_NW.t2p_interp];
    case 'BI'
        var2 = [W_NW.BI_logsmooth];
    case 'laminar'
        dist2input = [W_NW.inputDepths]; %input layer depths
        dist2input = dist2input(1:2:end);
        dist2input = dist2input - [W_NW.depth];% get distance to reversal - negative number is below reversal
        var2 = dist2input;
end

length(I)
var1 = var1(I);
var2 = var2(I);
dist2input = dist2input(I);


% I = find(BRI_NW<20);
%
% BRI_NW = BRI_NW(I);
%
% p2t_NW = [W_NW.t2p_interp];
% p2t_NW = p2t_NW(I);
%
% W_tmp = W_NW(I);
% figure(35); clf
%
% dist2input = [W_NW.inputDepths]; %input layer depths
% dist2input = dist2input(1:2:end);
% dist2input = dist2input - [W_NW.depth];% get distance to reversal - negative number is below reversal
% dist2input = dist2input(I);

figure(87); clf
scatter(var2, var1, 36, dist2input, 'filled'); hold on
set(gca,'yscale','log')
set(gca,'xscale','log')
colormap jet
colorbar

%set(gca, 'YScale', 'log')
%legend('Superficial Layer', 'Input Layer', 'Deep Layer')
title(tit)
xlabel(var2tit)
ylabel(var1tit)
hold off

%% Get W_keep - START HERE IF ALREADY GENERATED FILES

sess = [1:12 14:23 25 28:57];% 14:23 25 28:57];
dataPath = getpref('FREEVIEWING', 'PROCESSED_DATA_DIR');

meta_file = fullfile(fileparts(which('addFreeViewingPaths')), 'Data', 'datasets.csv');
% meta_file = fullfile(dataPath, 'datasets.xls');

data = readtable(meta_file);

W_keep2 = [];
for sessionId = sess
    
    if isnumeric(sessionId)
        thisSession = data(sessionId,:);
    end
    
    S.processedFileName = [thisSession.Tag{1} '.mat'];
    
    fname_new = fullfile(dataPath, 'waveform', 'waveform_new', S.processedFileName);
    
    load(fname_new);
    
    disp(['Getting spikes: session ' num2str(sessionId)])
    
    spfilename = [thisSession.Tag{1} '_kilo.mat'];
    spfname = fullfile(dataPath, 'spikes', spfilename);
    sp = load(spfname);
    sp = sp.sp; %CHANGE
    
    %     Keep1 = [W.isi_keep];
    %     Keep2 = [W.amp_keep];
    %     Keep = Keep1 + Keep2;
    tmpW = [];
    for g = 1:length(W)
        if contains(W(g).KSLabel, 'good')
            tmpW = [tmpW; W(g)];
        end
    end
    
    for i = 1:length(tmpW)
        tmpW(i).spikes = sp.st(sp.clu==tmpW(i).cid);
    end
    
    W_keep2 = [W_keep2; tmpW];
    
end

%% compute autocorr
savename = 'W_keep';
%W_keep = rmfield(W_keep,{'acor','ncor','pcor','tcor','normACOR'});

for i = 251:length(W_keep)
    
    [W_keep(i).tcor4,W_keep(i).acor4,W_keep(i).pcor4,W_keep(i).ncor4, W_keep(i).whichonesempty4] = comp_autocor_fast(W_keep(i).spikes,0.0004,0.040);
    [W_keep(i).tcor_smooth4,W_keep(i).acor_smooth4,W_keep(i).pcor_smooth4,W_keep(i).ncor_smooth4, ~,W_keep(i).BI,W_keep(i).PK,W_keep(i).LENO] = comp_autocor_fast_JM(W_keep(i).spikes,0.0004,0.040);
    [W_keep(i).tcor2,W_keep(i).acor2,W_keep(i).pcor2,W_keep(i).ncor2, ~] = comp_autocor_fast(W_keep(i).spikes,0.0002,0.040);
    [W_keep(i).tcor_smooth2,W_keep(i).acor_smooth2,W_keep(i).pcor_smooth2,W_keep(i).ncor_smooth2, ~,W_keep(i).BI2,W_keep(i).PK2,W_keep(i).LENO2] = comp_autocor_fast_JM(W_keep(i).spikes,0.0002,0.040);
    if mod(i, 50)==0
        save(savename, 'W_keep')
        disp(['Saved at index: ' num2str(i)])
    end
end

save(savename, 'W_keep')



