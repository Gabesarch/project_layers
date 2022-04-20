function eventTimes = getSaccadeEventTimes(Exp, which_noisetype)
% gets the event times of the current source density trials
% Inputs:
%   Exp              [struct] - Exp struct from io.dataFactoryGratingSubspace
% 
% jfm wrote it 2020
% ghs edit it 2020
% *** MODIFIED TO FLAG SACCADE ONSETS/OFFSETS INSTEAD OF FLASH ONSETS

if nargin < 2
  which_noisetype = 3;
end

sanitycheck = 0; % check spike histogram to verify CSD activity

Trials = length(Exp.D);
CSD = [];    % HERE THE STRUCT SAYS CSD, but time locking on SAC
CSD.Trials = [];
for k = 1:Trials
   ProtoName = Exp.D{k}.PR.name;
   %ProtoName
   type = 0;
   if (strcmp(ProtoName,'ForageProceduralNoise'))  % for Jake
       type = 1;
   end
   if (strcmp(ProtoName,'Forage'))   % for Shanna
       type = 2;
   end
  % type
   if (type > 0)
       NoiseType = Exp.D{k}.PR.noisetype;
       %disp(sprintf('Trial(%d); %s  NoiseType(%d)',k,ProtoName,NoiseType));
       if (type == 1) && (NoiseType == 3)
           CSD.Trials = [CSD.Trials ; [k 1]];  % 1 indicate contrast onset
       end
       if (type == 2) && ( (NoiseType == 3) || (NoiseType == 6) ) 
           if (NoiseType == 3) && (which_noisetype==3)
              CSD.Trials = [CSD.Trials ; [k 2]];
           end
           if (NoiseType == 6) && (which_noisetype==6)
              CSD.Trials = [CSD.Trials ; [k 3]];
           end
       end
   end
end

CSD.Onsets = [];
CSD.Types = [];
CSD.MoDir = [];
CSD.Offsets = [];
CSD.Onsets_Ephys = [];
CSD.Offsets_Ephys = [];
NTrials = length(CSD.Trials);
for k = 1:NTrials
    kk = CSD.Trials(k,1);
    type = CSD.Trials(k,2);
    
    %*** instead of NoiseHistory, step through Saccade History
    % NoHist = Exp.D{kk}.PR.NoiseHistory;
    slist = Exp.D{kk}.slist;            % saccade times relative to data
    eyesmo = Exp.D{kk}.eyeSmo;          % recorded eye position
    %********               
    
    %*** Noise History is time in column 1, and contrast (0 or 1) in col 2
    %** find all the Onsets, as transition from 0 to 1 in column 2
    for i = 1:size(slist,1)  % for i = 2:size(NoHist,1)
       if (1) % (NoHist(i-1,2) == 0) && (NoHist(i,2) >= 1)   % 
           tt = slist(i,1); % saccade onset time (float, in secs); % tt = NoHist(i,1);
           ott = slist(i,2); % saccade offset in secs
           CSD.Onsets = [CSD.Onsets ; tt];  % store Mat time
           CSD.Types = [CSD.Types ; type];  % 1 - contrast (Jake), 2 - contrast (Shanna), 3 - motion (Shanna)
           % CSD.MoDir = [CSD.MoDir ; NoHist(i,2)];
           % change MoDir to be the saccade direction
           if (1) % compute saccade direction here
             stasac = slist(i,4); % integer start time of saccade
             endsac = slist(i,5); % integer end time of saccade
             stx = eyesmo(stasac,2);
             sty = eyesmo(stasac,3);
             etx = eyesmo(endsac,2); 
             ety = eyesmo(endsac,3);
             sacvec = [etx-stx,ety-sty];
             ampo = norm(sacvec);
             ango = angle(complex(sacvec(1),sacvec(2))) * (360/(2*pi));
             if (ango < 0)
                 ango = ango + 360;
             end
             CSD.MoDir = ceil(ango/22.5);
           end
           %******* convert to Ephys time per trial start-clocks
           % NOTE, slist is time from trial start (while flashes were raw
           % matlab time, so subtraction step is not necessary here
           %**
           tt = tt + Exp.D{kk}.START_EPHYS;     % time from start of trial in ephys
           CSD.Onsets_Ephys = [CSD.Onsets_Ephys ; tt];
           %******
           CSD.Offsets = [CSD.Offsets ; ott];
           %******* convert to Ephys time per trial start-clocks
           ott = ott + Exp.D{kk}.START_EPHYS;     % time from start of trial in ephys
           CSD.Offsets_Ephys = [CSD.Offsets_Ephys ; ott];
       end 
    end
end

% Output (Other stuff was computed but we are ignoring that for now)
eventTimes = CSD.Onsets_Ephys;

if (sanitycheck == 1)
 for Unit = 1:size(Exp.sp,2) 
   NOnsets = length(CSD.Onsets_Ephys);
   SpkRaster = [];
   OffRaster = [];
   SpChan = Unit;
   for k = 1:NOnsets
     tt = CSD.Onsets_Ephys(k);
     % [min(Exp.sp{SpChan}.st),max(Exp.sp{SpChan}.st),tt]
     % input('check');
     z = find( (Exp.sp{SpChan}.st >= tt) & (Exp.sp{SpChan}.st < (tt+0.5) ) );
     if ~isempty(z)
      sptt = Exp.sp{SpChan}.st(z) - tt;  % time lock spikes relative to onset
      SpkRaster = [SpkRaster ; [(k*ones(size(sptt))) sptt]];
     end
     ott = CSD.Offsets_Ephys(k) - tt;
     OffRaster = [OffRaster ; [k ott]];
   end
   if ~isempty( SpkRaster )
     figure(10); hold off;
     h = plot(1000*OffRaster(:,2),OffRaster(:,1),'r.'); hold on;
     set(h,'Color',[1,0.8,0.8]);
     plot(1000*SpkRaster(:,2),SpkRaster(:,1),'k.'); hold on;
     xlabel('Time (ms)');
     ylabel('Trials');
     title(sprintf('CSD Raster: Unit(%d)',SpChan));
     input('check');
   else
     figure(10); hold off;
     
     disp(sprintf('Not plotting for unit %d, no spikes',Unit));
    end
 end
end
end

