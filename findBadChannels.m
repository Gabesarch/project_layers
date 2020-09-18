% filename1 = 'datasets.csv';
% % Read the CSV as a table
% t2 = readtable(filename1);
% % Add a new column to the end of the table
% numOfColumn = size(t, 2);
% newCol = repmat({'None'}, size(t,1), 1); %num2str(NaN(1, size(t,1)))'; % Your new column
% t.(numOfColumn) = newCol;
% % Change column name if needed
% t.Properties.VariableNames{numOfColumn} = 'deadChan';
% % Write to CSV file
% writetable(t, filename1)
%
% %%
%
% filename1 = 'datasets.csv';
% % Read the CSV as a table
% t = readtable(filename1);
% % Add a new column to the end of the table
% numOfColumn = size(t, 2);
% newCol = NaN(1, size(t,1))'; % Your new column
% t.(numOfColumn+1) = newCol;
% % Change column name if needed
% t.Properties.VariableNames{numOfColumn+1} = 'deadChan';
% % Write to CSV file
% writetable(t, filename1)

%%
% filename1 = 'datasets.csv';
% % Read the CSV as a table
% t2 = readtable(filename1);
% % Add a new column to the end of the table
% numOfColumn = size(t, 2);
% newCol = repmat({'None'}, size(t,1), 1); %num2str(NaN(1, size(t,1)))'; % Your new column
% t.(numOfColumn) = newCol;
% % Change column name if needed
% t.Properties.VariableNames{numOfColumn} = 'deadChan';
% % Write to CSV file
% writetable(t, filename1)
%
% %%
%
% filename1 = 'datasets.csv';
% % Read the CSV as a table
% t = readtable(filename1);
% % Add a new column to the end of the table
% numOfColumn = size(t, 2);
% newCol = NaN(1, size(t,1))'; % Your new column
% t.(numOfColumn+1) = newCol;
% % Change column name if needed
% t.Properties.VariableNames{numOfColumn+1} = 'deadChan';
% % Write to CSV file

%%
writetable(t, filename1)


%% Add Channels to existing based on bandpower
startOver = false;

filename1 = 'C:\Users\Gabe\Documents\MitchellLab\datasets.csv';
%filename1 = 'C:\Users\Gabe\Documents\MitchellLab\V1FreeViewingCode\Data\datasets.csv';

t = readtable(filename1);
numSess = size(t,1);

% left off at 36
for i = [28 32 37 38 45 46]%[30 33:36 39:42 44 47:48 50:51] %numSess
    disp('session: ')
    disp(i)
    csd.visualizeLayers(i)
    %[Exp, ~, lfp] = io.dataFactoryGratingSubspace(i);
    %p = csd.getLFPPower(lfp, [0 450], 'exclude', ~startOver, 'inds', [117000 118000]);
    
    x = 1;
    badChan = [];
    
    while ~isempty(x)
        prompt = 'Which values should you add? ';
        x = input(prompt);
        if ~isempty(x)
            badChan = [badChan x];
        end
    end
    badChan
    
    if isempty(badChan)
        %t(i, 29) = t(i, 29);
    else
        t(i, 29) = {num2str(sort(unique([str2num(t{i, 29}{:}) badChan])))};
    end
    
    
end


%% Find bad channels and add them to CSV file
filename1 = 'datasets.csv';
t = readtable(filename1);
numSess = size(t,1);

startOver = true;

for i = [28]% 32 37 38 45 46 49 52 53] %numSess
    disp('session: ')
    disp(i)
    [Exp, ~, lfp] = io.dataFactoryGratingSubspace(i);
    p = csd.getLFPPower(lfp, [0 450], 'exclude', ~startOver, 'inds', [117000 118000]);
    badChan = find(p < 0.2);
    badChan
    csd.plotLFP(lfp, [117000 118000], badChan)
    
    x = 1;
    if startOver == false
        badChan = [];
    end
    while ~isempty(x)
        if startOver
            prompt = 'Which values should you delete? ';
            x = input(prompt);
            if x == 100
                badChan = [];
                break
            end
            if ~isempty(x)
                badChan(badChan==x) = [];
            end
        elseif startOver == false
            prompt = 'Which values should you add? ';
            x = input(prompt);
            if ~isempty(x)
                badChan = [badChan x];
            end
        end
    end
    badChan
    
    if startOver
        if isempty(badChan)
            t(i, 29) = {'[]'};
        else
            t(i, 29) = {num2str(badChan)};
        end
    elseif startOver == false
        if isempty(badChan)
            %t(i, 29) = t(i, 29);
        else
            t(i, 29) = {num2str(sort(unique([str2num(t{i, 29}{:}) badChan])))};
        end
    end
    
    
end

%% Find bad channels and add them to CSV file
if ~exist('t')
    filename1 = 'datasets.csv';
    t = readtable(filename1);
end
numSess = size(t,1);

startOver = true;

for i = [ 52 53] %numSess
    disp('session: ')
    disp(i)
    [Exp, ~, lfp] = io.dataFactoryGratingSubspace(i);
    p = csd.getLFPPower(lfp, [0 450], 'exclude', ~startOver, 'inds', [117000 118000]);
    badChan = [];
    badChan
    csd.plotLFP(lfp, [117000 118000], badChan)
    
    x = 1;
    while ~isempty(x)
        prompt = 'Which values should you add? ';
        x = input(prompt);
        if ~isempty(x)
            badChan = [badChan x];
        end
    end
    
    badChan
    
    t(i, 29) = {num2str(sort(unique([str2num(t{i, 29}{:}) badChan])))};
    
    
    
end