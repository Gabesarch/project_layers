function varargout = addMarmoPipe()
% addMarmoPipe adds the required paths for MarmoPipe. It assumes that all
% paths are relative to this function

% where is marmoPipe?
marmoPath = fileparts(mfilename('fullpath'));


% dirs{n}='Import'; % MarmoPipe Import tools

% Add the external repositories we need
directory = fullfile(marmoPath, '');
dirs = {'MarmosetScripts', ...
    'ephys-matlab', ...
    '+csd'};

% remove directories that could not be found
remove = cellfun(@(x) ~exist(fullfile(directory, x), 'dir'), dirs);

if any(remove) % throw warning
    removeList = dirs(remove);
    dirs(remove) = [];
    fprintf('Could not find %d external repositories. \nI''m trusting you to add them.\nListing here:\n', numel(removeList))
    
    for i = 1:numel(removeList)
        fprintf('%s\n', removeList{i});
    end
end
    
addPathsWithoutGit(directory, dirs);

% add MarmoPipe and Import
%addPathsWithoutGit(marmoPath, {'Import', 'MarmoPipe'})

if nargout > 0 % if argument is requested, pass out the path to marmoPath
    varargout{1} = marmoPath;
end

function addPathsWithoutGit(directory, dirList)
% addPathsWithoutGit is essentially the same as addpath(genpath(dirList)),
% but it ignores the .git hidden folders

for j=1:length(dirList)
    a=genpath(fullfile(directory, dirList{j}));
    % convert list of paths into a cellarray of subfolders
    b=textscan(a,'%s','delimiter',';');
    b=b{1};
    % remove git hidden folders
    b(~cellfun(@isempty,strfind(b,'.git')))=[]; %#ok<STRCLFH>
    % add all the paths
    addpath(b{:})
    fprintf('added [%s] to the path\n', dirList{j});
end

cd('C:\Users\Gabe\Documents\MitchellLab\V1FreeViewingCode')
addpath('Data')
addFreeViewingPaths('gabelaptop')