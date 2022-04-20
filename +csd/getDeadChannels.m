function lfp = getDeadChannels(lfp,Exp,exclude_threshold)



unique_x = unique(Exp.osp.xcoords);
num_shanks = size(unique_x, 1);
if nargin < 3
    exclude_threshold = 5*num_shanks; % if more than this just leave alone
end
ycoords = Exp.osp.ycoords;
diff_ycoords = abs(mode(diff(Exp.osp.ycoords)));
shank_len = size(lfp.data, 2)/num_shanks;
ycoords_new = repmat(flip(linspace(diff_ycoords, diff_ycoords*shank_len, shank_len)), 1, 1)';
if num_shanks>1
    [p, max_diff] = max(diff(ycoords));
    ycoords1 = ycoords(1:max_diff);
    exclude = [];
    for i = 1:length(ycoords_new)
        y_cur = ycoords_new(i);
        contains = any(ycoords1(:)==y_cur);
        if ~contains
            exclude = [exclude i];
        end
    end
    ycoords2 = ycoords(max_diff+1:end);
    for i = 1:length(ycoords_new)
        y_cur = ycoords_new(i);
        contains = any(ycoords2(:)==y_cur);
        if ~contains
            exclude = [exclude i+32];
        end
    end
else
    [p, max_diff] = max(diff(ycoords));
    ycoords1 = ycoords(1:max_diff);
    exclude = [];
    for i = 1:length(ycoords_new)
        y_cur = ycoords_new(i);
        contains = any(ycoords1(:)==y_cur);
        if ~contains
            exclude = [exclude i];
        end
    end
end

disp('Number of dead channels identified:')
disp(length(exclude))
if length(exclude)<exclude_threshold
    lfp.deadChan = exclude;
else
    warning('Too many dead channels exist. Leaving them in to avoid too much removal')
    lfp.deadChan = [];
end


end

