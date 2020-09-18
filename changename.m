

sess = 1:57 %1:57;
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
    fname = fullfile(dataPath, 'spikes', S.processedFileName);
    
    fname_new = [thisSession.Tag{1} '_kilo.mat'];
    
    movefile(S.processedFileName,fname_new);
    
    catch
        error = [error sessionId];
    end
end