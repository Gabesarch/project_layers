function [tcor,acor,pcor,ncor] = comp_autocor_fast_with_intervals(sptimes,binsize,maxlag,intervals)
% function [auto] = comp_autocor_fast(sptimes,binsize,maxlag)
%
%** inputs:  sptimes, Nx1 list of spike times in secs
%                if empty, it will build an artificial spike train
%                    with a refractory and burst process
%            binsize, 0.2 ms default, size of bins for spikes
%            maxlag, maximum lag to show autocor out to, 30ms default
%            intervals:   a Mx2 array, where M is the number of intervals
%                         and the first column is a start time (in secs)
%                         and the second column is a end time (in secs)
%                  ... it will only analyze spikes inside the intervals
%***
%*** outputs:  returns autocorrelation, tcor is timelags (0 to maxlag), 
%***                                    acor is autocorrelation 
%***                                    pcor is poisson expectation of auto
%***                                    ncor (acor ./ pcor)

  if isempty(sptimes)
      disp('Building example spike train ...');
      GoTime = 0.0;
      FiTime = 100.0;   % simulate 1000 seconds of artifical process
      %****** First, build poisson spike from random drawn times
      Rate = 100;         % mean rate for a Poisson process (random spikes)
      NSpikes = (Rate * FiTime);
      sptimes = GoTime + (FiTime - GoTime) * rand(NSpikes,1);
      sptimes = sort(sptimes);
      %****** Apply a 1.5 ms refractory period to process
      xtimes = [1 ; abs(diff(sptimes))];
      zz = find( xtimes < 0.0015);
      sptimes(zz) = sptimes(zz) + 0.0015;  
      % shift refractory collisions 1ms into the future, which
      % should make a burst like bump at 1.5 - 3 ms
      disp('... finished example train');
  end
                                          
  if isempty(binsize)
      BinSize = 0.0004;  % 0.4 ms (in secs)
  else
      BinSize = binsize;
  end 
  if isempty(maxlag)
      MaxLag = 0.040;    % 40 ms lag (in secs)
  else
      MaxLag = maxlag;
  end
  NLag = 1+floor(MaxLag/BinSize);
  tcor = (1000 * BinSize) * (0:(NLag-1));
  
  %**** reformulate sptimes so it only includes intervals
  if ~isempty(intervals)
     newsptimes = [];
     for k = 1:size(intervals,1)
         zz = find( (sptimes >= intervals(k,1)) & (sptimes < intervals(k,2)) );
         newsptimes = [newsptimes ; sptimes(zz)];  
     end
     if (isempty(newsptimes))
         disp('Failure to identify any spikes within intervals specified');
         acor = nan(1,NLag);
         pcor = acor;
         ncor = acor;
         return;
     else
        sptimes = newsptimes; 
     end
  end
  
  %******* use xcorr to compute autocorrelation on binned spikes
  disp('Computing correlation over lags');
  Asum = zeros(1,NLag);
  k10 = floor(length(sptimes)/10);  % report progress in 10ths
  for k = 1:length(sptimes)
      spt = sptimes(k);
      %** any time zero multiplies, you get zero
      %** thus you only need those moments that bin 0 has a one
      zz = find( (sptimes > spt) & (sptimes < (spt+MaxLag+BinSize)) );
      if ~isempty(zz)
          tlags = 1 + floor( (sptimes(zz)-spt)/BinSize );
          Asum(tlags) = Asum(tlags)+1;
      end
      %***  report progress on command line
      if (mod(k,k10) == 0)
          disp(sprintf('Progress %d percent',ceil(k*100/length(sptimes))));
      end
      %*******
  end
  acor = Asum / length(sptimes);
  %****** normalize by expectation for a Poisson process (flat autocor)
  disp('Computing Poisson expectation');
  if ~isempty(intervals)
      totime = 0;
      for k = 1:length(intervals)
          totime = totime + (intervals(k,2) - intervals(k,1));
      end
  else
      totime = max(sptimes) - min(sptimes);
  end
  rate = length(sptimes)/totime;
  arate = (rate*binsize);
  pcor = arate * ones(size(acor));
  ncor = acor ./ pcor;
  
  %******* if you want to plot the result
  if (1)
    figure;
    plot(tcor,acor,'k.-'); hold on;
    plot(tcor,pcor,'r-');
    xlabel('Time (ms)');
    ylabel('Autocor');
  end
  
return;