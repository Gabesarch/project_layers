
error = [];
for i = 1:57
    try 
        %stats = csd.dataFactoryCSD(i, 'type', 'csd');
        stats2 = csd.dataFactoryCSD(i, 'type', 'gamma');%, 'method', 'weightedMin', 'plotIt', false);
        %gamma = csd.getGamma(lfp, 'plotIt', true, 'method', 'weightedMin')
    catch
        error = [error i];
    end
end