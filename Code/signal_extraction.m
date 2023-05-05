%Elongation Signal Extraction

%signal_extraction(hbmlp13,noise offset,gradient)

%DV for DV gradient
%AP for AP gradient
function [hbmlp13, alert] = signal_extraction(hbmlp13,insert,varargin)
%use for control
if isempty(varargin)
    noNoise = 2;
    distance = 0;
    gradient = 'AP';
    globalthresh = 0;
elseif size(varargin,2) == 1
    noNoise = varargin{1};
    distance = 0;
    gradient = 'AP';
    globalthresh = 0;
elseif size(varargin,2) == 2
    noNoise = varargin{1};
    gradient = varargin{2};
    distance = 0;
    globalthresh = 0;
elseif size(varargin,2) == 3
    noNoise = varargin{1};
    gradient = varargin{2};
    globalthresh = varargin{3};
    distance = 0;
else
    noNoise = varargin{1};
    gradient = varargin{2};
    distance = varargin{4};
    globalthresh = varargin{3};

end

%insert 
%mlp = 4.5kb
%mglp = 6kb
%mglGp = 8.5kb


%% Adjust Signal, find active
for j = 1:length(hbmlp13)
    hbmlp13(j).Mm = double(hbmlp13(j).Mm); %convert to double
    hbmlp13(j).Pm = double(hbmlp13(j).Pm);
    for i = 1:length(hbmlp13(j).M)     
        hbmlp13(j).M(:,i) = hbmlp13(j).M(:,i)-hbmlp13(j).Mm(:,i); %subtract mean value of nuclei first
        hbmlp13(j).M(:,i) = hbmlp13(j).M(:,i)-min(hbmlp13(j).M(:,i));
        hbmlp13(j).P(:,i) = hbmlp13(j).P(:,i)-hbmlp13(j).Pm(:,i); %subtract mean value of nuclei first
        hbmlp13(j).P(:,i) = hbmlp13(j).P(:,i)-min(hbmlp13(j).P(:,i));
    end
end

for j = 1:length(hbmlp13)

    [hbmlp13(j).intT] = linspace(hbmlp13(j).T(1),hbmlp13(j).T(end));
    [hbmlp13(j).trueM] = interp1(hbmlp13(j).T,hbmlp13(j).M,hbmlp13(j).intT);
    [hbmlp13(j).trueP] = interp1(hbmlp13(j).T,hbmlp13(j).P,hbmlp13(j).intT);
    
    for i = 1:length(hbmlp13(j).M)     
        
        %try smoothing
        hbmlp13(j).M(:,i) = smooth(hbmlp13(j).M(:,i),4);
        hbmlp13(j).P(:,i) = smooth(hbmlp13(j).P(:,i),4);

        
        sortM{j}(:,i) = sort(hbmlp13(j).M(:,i));
        sortP{j}(:,i) = sort(hbmlp13(j).P(:,i));
        maxM{j} = sortM{j}(end-2,:);
        maxP{j} = sortP{j}(end-2,:);
        threshM(j) = .8*mean(mean(hbmlp13(j).Mm));
        threshP(j) = .8*mean(mean(hbmlp13(j).Pm));
       
    end
    [hbmlp13(j).threshM] = threshM(j);
    [hbmlp13(j).threshP] = threshP(j);
    
    [hbmlp13(j).active] = find(maxM{j} > threshM(j) & maxP{j} > threshP(j));
    [hbmlp13(j).intM] = interp1(hbmlp13(j).T,hbmlp13(j).M,hbmlp13(j).intT);
    [hbmlp13(j).intP] = interp1(hbmlp13(j).T,hbmlp13(j).P,hbmlp13(j).intT);
    
    for i = 1:length(hbmlp13(j).M)
        
    [hbmlp13(j).normM(:,i)] = hbmlp13(j).trueM(:,i)./max(hbmlp13(j).trueM(:,i));
    [hbmlp13(j).normP(:,i)] = hbmlp13(j).trueP(:,i)./max(hbmlp13(j).trueP(:,i));
    
    [hbmlp13(j).normintM(:,i)] = hbmlp13(j).intM(:,i)./max(hbmlp13(j).intM(:,i));
    [hbmlp13(j).normintP(:,i)] = hbmlp13(j).intP(:,i)./max(hbmlp13(j).intP(:,i));
    end
    
    [hbmlp13(j).activeM] = find(maxM{j} > threshM(j));
    if isempty(hbmlp13(j).activeM)
        [hbmlp13(j).activeP] = [];
    else
        for i = 1:length(hbmlp13(j).activeM)
            actPtmp(i) = maxP{j}(i) > threshP(j);
        end
        [hbmlp13(j).activeP] = find(actPtmp); %values are indices of hbmlp13(j).activeM not of all nuclei. Need to refer to .activeM to find actual nucleus 
    end
    
%     [hbmlp13(j).inactive] = find(maxM{j} < threshM(j) & maxP{j} > threshP(j));
    %find how many nuclei become active over time
    if isempty(hbmlp13(j).active)
        [hbmlp13(j).activetime] = [];
    else
        for i = 1:length(hbmlp13(j).active)
            T_M = find(hbmlp13(j).intM(:,hbmlp13(j).active(i)) > hbmlp13(j).threshM); % find when signal is "active"
            T_active(i) = hbmlp13(j).intT(T_M(2)); %find second interpolated time point where signal is above active threshold
        end
        nuc = 0;
        for i = 1:length(hbmlp13(j).intT)
            nuc_active = find(T_active == hbmlp13(j).intT(i));
            if i ==1
                nuc(i) = nuc + size(nuc_active,2);
            else 
                nuc(i) = nuc(i-1) +size(nuc_active,2);
            end
        end
        [hbmlp13(j).activetime] = nuc;
        [hbmlp13(j).activenuc] = T_active;
        clear T_M T_active actPtmp
    end
    [hbmlp13(j).normT] = hbmlp13(j).intT./hbmlp13(j).intT(end);
end

%% Below looks at the distance between max pixels of MS2 and PP7
% if they are more than 3 pixels away in a given time point it is flagged
% if there are 2 time points with flags the nuclei is no longer considered
% active

if distance == 1
for i = 1:length(hbmlp13)
    if isempty(hbmlp13(i).active)
        continue
    else
        for j = 1:length(hbmlp13(i).active)
            for q = 1:length(hbmlp13(i).T)
                [hbmlp13(i).dist{j}(q,1)] = pdist([hbmlp13(i).MS2pix(q,:,hbmlp13(i).active(j)); hbmlp13(i).PP7pix(q,:,hbmlp13(i).active(j))]);
                if hbmlp13(i).M(q,hbmlp13(i).active(j)) >= threshM & hbmlp13(i).P(q,hbmlp13(i).active(j)) >= threshP
                    if hbmlp13(i).dist{j}(q,1) < 5
                        [hbmlp13(i).dist{j}(q,2)] = 0;
                    else
                        [hbmlp13(i).dist{j}(q,2)] = 1;
                        
                    end
                else
                    [hbmlp13(i).dist{j}(q,2)] = NaN;
                end
            end
            if nansum(hbmlp13(i).dist{j}(:,2) >= 1)
                alert{i} = hbmlp13(i).active(j);
            end
        end

    end
end
for i = 1:length(hbmlp13)
    if isempty(alert{i})
        continue
    else
        for q = 1:length(alert{i})
            var(q) = intersect(hbmlp13(i).active,alert{i});
            hbmlp13(i).active(hbmlp13(i).active == var(q)) = [];
        end
    end
    clear var
end
end

%% Transctiption Time
%{
if globalthresh == 0
A = .3; %Percent of max to determine elongation time
Pdelay = noNoise; %Time steps of delay to remove noise in beginning of cycle
for j = 1:length(hbmlp13)
    for i = 1:length(hbmlp13(j).active)
        %Find when MS2 is active first
        tmpM = find(hbmlp13(j).intM(:,hbmlp13(j).active(i)) > A*max(hbmlp13(j).M(Pdelay:end,hbmlp13(j).active(i)))); %find when MS2 is above threshold
        if tmpM(1) == 1
            if hbmlp13(j).intM(1:Pdelay(1),hbmlp13(j).active(i)) > A*max(hbmlp13(j).M(Pdelay:end,hbmlp13(j).active(i))) %If signal is high from beginning
                intM = hbmlp13(j).intT(1);
            else
                tmpM = find(hbmlp13(j).intM(Pdelay(1):-1:1,hbmlp13(j).active(i)) < A*max(hbmlp13(j).M(Pdelay:end,hbmlp13(j).active(i)))); %If signal is already on step backwards until you find where it crosses the threshold
                intM = interp1([hbmlp13(j).intM(Pdelay(1)-tmpM(1),hbmlp13(j).active(i)) hbmlp13(j).intM(Pdelay(1)-tmpM(1)+1,hbmlp13(j).active(i))],[hbmlp13(j).intT(Pdelay(1)-tmpM(1)) hbmlp13(j).intT(Pdelay(1)-tmpM(1)+1)],A*max(hbmlp13(j).M(Pdelay:end,hbmlp13(j).active(i))));
            end
        else
            intM = interp1([hbmlp13(j).intM(tmpM(1)-1,hbmlp13(j).active(i)) hbmlp13(j).intM(tmpM(1),hbmlp13(j).active(i))],[hbmlp13(j).intT(tmpM(1)-1) hbmlp13(j).intT(tmpM(1))],A*max(hbmlp13(j).M(Pdelay:end,hbmlp13(j).active(i))));
        end
        
        %now look for PP7 after MS2 is active
        Pstart = find(hbmlp13(j).intT > intM);
        tmpP = find(hbmlp13(j).intP(Pstart(1):end,hbmlp13(j).active(i)) > A*max(hbmlp13(j).P(Pdelay:end,hbmlp13(j).active(i)))); %look at signal after MS2 is active
       
        
        if isempty(tmpP)
            intP = NaN;
        elseif tmpP(1) == 1
            tmpP = find(hbmlp13(j).intP(Pstart(1):-1:1,hbmlp13(j).active(i)) < A*max(hbmlp13(j).P(Pdelay:end,hbmlp13(j).active(i)))); %walk backwards to find when signal became active
            if isempty(tmpP)
                intP = NaN;
            else
%                 intP = interp1([hbmlp13(j).intP(Pstart(1)-tmpP(1)+1,hbmlp13(j).active(i)) hbmlp13(j).intP(Pstart(1)-tmpP(1)+2,hbmlp13(j).active(i))],[hbmlp13(j).intT(Pstart(1)-tmpP(1)+1) hbmlp13(j).intT(Pstart(1)-tmpP(1)+2)],A*max(hbmlp13(j).P(Pdelay:end,hbmlp13(j).active(i))));
                intP = NaN; %ignore nuclei that give negative values
            end
        else
            intP = interp1([hbmlp13(j).intP(Pstart(1)+tmpP(1)-2,hbmlp13(j).active(i)) hbmlp13(j).intP(Pstart(1)+tmpP(1)-1,hbmlp13(j).active(i))], [hbmlp13(j).intT(tmpP(1)+Pstart(1)-2) hbmlp13(j).intT(tmpP(1)+Pstart(1)-1)],A*max(hbmlp13(j).P(Pdelay:end,hbmlp13(j).active(i))));

        end
        
        [hbmlp13(j).actM(i)] = intM;
        [hbmlp13(j).actP(i)] = intP;
        [hbmlp13(j).dT(i)] = intP - intM;
        [hbmlp13(j).rate(i)] = insert./hbmlp13(j).dT(i);
    end
end
end
%}

%% Transctiption Time Normalized


if globalthresh == 0
A = .3; %Percent of max to determine elongation time
Pdelay = noNoise; %Time steps of delay to remove noise in beginning of cycle
for j = 1:length(hbmlp13)
    for i = 1:length(hbmlp13(j).active)
        %Find when MS2 is active first
        tmpM = find(hbmlp13(j).normintM(:,hbmlp13(j).active(i)) > A*max(hbmlp13(j).normintM(Pdelay:end,hbmlp13(j).active(i)))); %find when MS2 is above threshold
        if tmpM(1) == 1
            if hbmlp13(j).normintM(1:Pdelay(1),hbmlp13(j).active(i)) > A*max(hbmlp13(j).normintM(Pdelay:end,hbmlp13(j).active(i))) %If signal is high from beginning
                normintM = hbmlp13(j).intT(1);
            else
                tmpM = find(hbmlp13(j).normintM(Pdelay(1):-1:1,hbmlp13(j).active(i)) < A*max(hbmlp13(j).normintM(Pdelay:end,hbmlp13(j).active(i)))); %If signal is already on step backwards until you find where it crosses the threshold
                normintM = interp1([hbmlp13(j).normintM(Pdelay(1)-tmpM(1),hbmlp13(j).active(i)) hbmlp13(j).normintM(Pdelay(1)-tmpM(1)+1,hbmlp13(j).active(i))],[hbmlp13(j).intT(Pdelay(1)-tmpM(1)) hbmlp13(j).intT(Pdelay(1)-tmpM(1)+1)],A*max(hbmlp13(j).normintM(Pdelay:end,hbmlp13(j).active(i))));
            end
        else
            normintM = interp1([hbmlp13(j).normintM(tmpM(1)-1,hbmlp13(j).active(i)) hbmlp13(j).normintM(tmpM(1),hbmlp13(j).active(i))],[hbmlp13(j).intT(tmpM(1)-1) hbmlp13(j).intT(tmpM(1))],A*max(hbmlp13(j).normintM(Pdelay:end,hbmlp13(j).active(i))));
        end
        
        %now look for PP7 after MS2 is active
        Pstart = find(hbmlp13(j).intT > normintM);
        tmpP = find(hbmlp13(j).normintP(Pstart(1):end,hbmlp13(j).active(i)) > A*max(hbmlp13(j).normintP(Pdelay:end,hbmlp13(j).active(i)))); %look at signal after MS2 is active
       
        
        if isempty(tmpP)
            normintP = NaN;
        elseif tmpP(1) == 1
            tmpP = find(hbmlp13(j).normintP(Pstart(1):-1:1,hbmlp13(j).active(i)) < A*max(hbmlp13(j).normintP(Pdelay:end,hbmlp13(j).active(i)))); %walk backwards to find when signal became active
            if isempty(tmpP)
                normintP = NaN;
            else
%                 normintP = interp1([hbmlp13(j).normintP(Pstart(1)-tmpP(1)+1,hbmlp13(j).active(i)) hbmlp13(j).normintP(Pstart(1)-tmpP(1)+2,hbmlp13(j).active(i))],[hbmlp13(j).intT(Pstart(1)-tmpP(1)+1) hbmlp13(j).intT(Pstart(1)-tmpP(1)+2)],A*max(hbmlp13(j).P(Pdelay:end,hbmlp13(j).active(i))));
                normintP = NaN; %ignore nuclei that give negative values
            end
        else
            normintP = interp1([hbmlp13(j).normintP(Pstart(1)+tmpP(1)-2,hbmlp13(j).active(i)) hbmlp13(j).normintP(Pstart(1)+tmpP(1)-1,hbmlp13(j).active(i))], [hbmlp13(j).intT(tmpP(1)+Pstart(1)-2) hbmlp13(j).intT(tmpP(1)+Pstart(1)-1)],A*max(hbmlp13(j).normintP(Pdelay:end,hbmlp13(j).active(i))));

        end
        
        [hbmlp13(j).actM(i)] = normintM;
        [hbmlp13(j).actP(i)] = normintP;
        [hbmlp13(j).dT(i)] = normintP - normintM;
        [hbmlp13(j).rate(i)] = insert./hbmlp13(j).dT(i);
    end
end
end


%% Global threshold based off % of average max
if globalthresh == 1
A = 0.3; %percent of max
Pdelay = noNoise; %Time steps of delay to remove noise in beginning of cycle
for j = 1:length(hbmlp13)
    MS2max(j) = A*mean(sortM{j}(end,hbmlp13(j).active),2); %get point to measure time difference. take percent of average max signal of active nuclei
    PP7max(j) = A*mean(sortP{j}(end,hbmlp13(j).active),2);
    
    [hbmlp13(j).active] = find(maxM{j} > MS2max(j) & maxP{j} > PP7max(j)); %new threshold when using global timing threshold (signal needs to be above the threshold)
    for i = 1:length(hbmlp13(j).active)
        %Find when MS2 is active first
        tmpM = find(hbmlp13(j).normintM(:,hbmlp13(j).active(i)) > MS2max(j)); %find when MS2 is above threshold
        if tmpM(1) == 1
            if hbmlp13(j).normintM(1:Pdelay(1),hbmlp13(j).active(i)) > MS2max(j) %If signal is high from beginning
                normintM = hbmlp13(j).intT(1);
            else
                tmpM = find(hbmlp13(j).normintM(Pdelay(1):-1:1,hbmlp13(j).active(i)) < MS2max(j)); %If signal is already on step backwards until you find where it crosses the threshold
                normintM = interp1([hbmlp13(j).normintM(Pdelay(1)-tmpM(1),hbmlp13(j).active(i)) hbmlp13(j).normintM(Pdelay(1)-tmpM(1)+1,hbmlp13(j).active(i))],[hbmlp13(j).intT(Pdelay(1)-tmpM(1)) hbmlp13(j).intT(Pdelay(1)-tmpM(1)+1)],MS2max(j));
            end
        else
            normintM = interp1([hbmlp13(j).normintM(tmpM(1)-1,hbmlp13(j).active(i)) hbmlp13(j).normintM(tmpM(1),hbmlp13(j).active(i))],[hbmlp13(j).intT(tmpM(1)-1) hbmlp13(j).intT(tmpM(1))],MS2max(j));
        end
        
        %now look for PP7 after MS2 is active
        Pstart = find(hbmlp13(j).intT > normintM);
        tmpP = find(hbmlp13(j).normintP(Pstart(1):end,hbmlp13(j).active(i)) > PP7max(j)); %look at signal after MS2 is active
       
        
        if isempty(tmpP)
            normintP = NaN;
        elseif tmpP(1) == 1
            tmpP = find(hbmlp13(j).normintP(Pstart(1):-1:1,hbmlp13(j).active(i)) < PP7max(j)); %walk backwards to find when signal became active
            if isempty(tmpP)
                normintP = NaN;
            else
%                 normintP = interp1([hbmlp13(j).normintP(Pstart(1)-tmpP(1)+1,hbmlp13(j).active(i)) hbmlp13(j).normintP(Pstart(1)-tmpP(1)+2,hbmlp13(j).active(i))],[hbmlp13(j).intT(Pstart(1)-tmpP(1)+1) hbmlp13(j).intT(Pstart(1)-tmpP(1)+2)],A*max(hbmlp13(j).P(Pdelay:end,hbmlp13(j).active(i))));
                normintP = NaN; %ignore nuclei that give negative values
            end
        else
            normintP = interp1([hbmlp13(j).normintP(Pstart(1)+tmpP(1)-2,hbmlp13(j).active(i)) hbmlp13(j).normintP(Pstart(1)+tmpP(1)-1,hbmlp13(j).active(i))], [hbmlp13(j).intT(tmpP(1)+Pstart(1)-2) hbmlp13(j).intT(tmpP(1)+Pstart(1)-1)],PP7max(j));

        end
        
        [hbmlp13(j).actM(i)] = normintM;
        [hbmlp13(j).actP(i)] = normintP;
        [hbmlp13(j).dT(i)] = normintP - normintM;
        [hbmlp13(j).rate(i)] = insert./hbmlp13(j).dT(i);
        
    end
    [hbmlp13(j).MS2max(j)] = MS2max(j);
    [hbmlp13(j).PP7max(j)] = PP7max(j);
end
end
%%
%{
%% Transcription Time
A = .3; %Percent of max to determine elongation time
Pdelay = noNoise; %Time steps of delay to remove noise in beginning of cycle
for j = 1:length(hbmlp13)
    min5 = find(hbmlp13(j).intT > noNoise); %find for all replicates index when intT = 2. Use to look after noise for P signal

    for i = 1:length(hbmlp13(j).active)
        tmpM = find(hbmlp13(j).intM(:,hbmlp13(j).active(i)) > A*max(hbmlp13(j).M(Pdelay:end,hbmlp13(j).active(i))));    %take max of interp signal or real signal?????!!!!!!!
        tmpP = find(hbmlp13(j).intP(min5(1):end,hbmlp13(j).active(i)) > A*max(hbmlp13(j).P(Pdelay:end,hbmlp13(j).active(i)))); %look at signal after 4.6 min
        
        %If signal is already above 30% after 2 min
        if isempty(tmpP) || tmpP(1) == 1 
            if hbmlp13(j).intP(1:min5(1),hbmlp13(j).active(i)) > A*max(hbmlp13(j).P(Pdelay:end,hbmlp13(j).active(i))) %If signal is high from beginning
                intP = hbmlp13(j).intT(1);
            else
            tmpP = find(hbmlp13(j).intP(min5(1):-1:1,hbmlp13(j).active(i)) < A*max(hbmlp13(j).P(Pdelay:end,hbmlp13(j).active(i)))); %walk backwards to find when signal became active
            intP = interp1([hbmlp13(j).intP(min5(1)-tmpP(1)+1,hbmlp13(j).active(i)) hbmlp13(j).intP(min5(1)-tmpP(1)+2,hbmlp13(j).active(i))],[hbmlp13(j).intT(min5(1)-tmpP(1)+1) hbmlp13(j).intT(min5(1)-tmpP(1)+2)],A*max(hbmlp13(j).P(Pdelay:end,hbmlp13(j).active(i))));
            end
        
%         %Disregard Single spikes and skip to next peak above 30%
%         elseif length(tmpP) > 6 && tmpP(5) ~= tmpP(1)+4
%             tmpP = find(hbmlp13(j).intP(1:tmpP(5)+min5(1)-1,hbmlp13(j).active(i)) < A*max(hbmlp13(j).P(Pdelay:end,hbmlp13(j).active(i)))); %find point just before actual activation           
%             intP = interp1([hbmlp13(j).intP(tmpP(end),hbmlp13(j).active(i)) hbmlp13(j).intP(tmpP(end)+1,hbmlp13(j).active(i))],[hbmlp13(j).intT(tmpP(end)) hbmlp13(j).intT(tmpP(end)+1)],A*max(hbmlp13(j).P(Pdelay:end,hbmlp13(j).active(i))));
%             
        else
            intP = interp1([hbmlp13(j).intP(tmpP(1)+min5(1)-2,hbmlp13(j).active(i)) hbmlp13(j).intP(tmpP(1)+min5(1)-1,hbmlp13(j).active(i))],[hbmlp13(j).intT(tmpP(1)+min5(1)-2) hbmlp13(j).intT(tmpP(1)+min5(1)-1)],A*max(hbmlp13(j).P(Pdelay:end,hbmlp13(j).active(i))));
        end
        
        if tmpM(1) == 1
            if hbmlp13(j).intM(1:min5(1),hbmlp13(j).active(i)) > A*max(hbmlp13(j).M(Pdelay:end,hbmlp13(j).active(i))) %If signal is high from beginning
                intM = hbmlp13(j).intT(1);
            else
                tmpM = find(hbmlp13(j).intM(min5(1):-1:1,hbmlp13(j).active(i)) < A*max(hbmlp13(j).M(Pdelay:end,hbmlp13(j).active(i))));
                intM = interp1([hbmlp13(j).intM(min5(1)-tmpM(1)+1,hbmlp13(j).active(i)) hbmlp13(j).intM(min5(1)-tmpM(1)+2,hbmlp13(j).active(i))],[hbmlp13(j).intT(min5(1)-tmpM(1)+1) hbmlp13(j).intT(min5(1)-tmpM(1)+2)],A*max(hbmlp13(j).M(Pdelay:end,hbmlp13(j).active(i))));
            end
        else
            intM = interp1([hbmlp13(j).intM(tmpM(1)-1,hbmlp13(j).active(i)) hbmlp13(j).intM(tmpM(1),hbmlp13(j).active(i))],[hbmlp13(j).intT(tmpM(1)-1) hbmlp13(j).intT(tmpM(1))],A*max(hbmlp13(j).M(Pdelay:end,hbmlp13(j).active(i))));
        end
        
        [hbmlp13(j).actM(i)] = intM;
        [hbmlp13(j).actP(i)] = intP;
        [hbmlp13(j).dT(i)] = intP - intM;
%         [hbmlp13(j).dT(i)] = hbmlp13(j).intT(tmpP(1)+35) - hbmlp13(j).intT(tmpM(1));
        
    end
end
%}
%% Bin nuclei
for j = 1:length(hbmlp13)
    timeframe = round(length(hbmlp13(j).T)/2);
    bin_lines = [1:32:512 512];
    for i = 1:16
        if strcmp(gradient,'AP') == 1
            [hbmlp13(j).bin{i}] = find(hbmlp13(j).lineage_cx(timeframe,:) > bin_lines(i) & hbmlp13(j).lineage_cx(timeframe,:) < bin_lines(i+1));
        elseif strcmp(gradient,'DV') == 1
            [hbmlp13(j).bin{i}] = find(hbmlp13(j).lineage_cy(timeframe,:) > bin_lines(i) & hbmlp13(j).lineage_cy(timeframe,:) < bin_lines(i+1));
        end
    end
end

%% Output (using interpolated trace)

%Ouput based off bin
for j = 1:length(hbmlp13)
    dX(j) = hbmlp13(j).intT(2)-hbmlp13(j).intT(1);
    for i = 1:16
        if isempty(hbmlp13(j).bin{i})
            tmpOutM{i} = 0;
            tmpOutP{i} = 0;
            
        else
            for q = 1:length(hbmlp13(j).bin{i})
               tmpOutM{i}(q) = trapz(hbmlp13(j).intT,hbmlp13(j).normintM(:,hbmlp13(j).bin{i}(q)))/(dX(j)*100);
               tmpOutP{i}(q) = trapz(hbmlp13(j).intT,hbmlp13(j).normintP(:,hbmlp13(j).bin{i}(q)))/(dX(j)*100);
            end
        end
        [hbmlp13(j).OutputBinM{i}] = tmpOutM{i};
        [hbmlp13(j).OutputBinP{i}] = tmpOutP{i};
    end
end

%Output
for j = 1:length(hbmlp13)
    dX(j) = hbmlp13(j).intT(2)-hbmlp13(j).intT(1);
    for i = 1:length(hbmlp13(j).active)
        [hbmlp13(j).OutputM(i)] = trapz(hbmlp13(j).intT,hbmlp13(j).normintM(:,hbmlp13(j).active(i)))/(dX(j)*100);
        [hbmlp13(j).OutputP(i)] = trapz(hbmlp13(j).intT,hbmlp13(j).normintP(:,hbmlp13(j).active(i)))/(dX(j)*100);
    end
end


%Look at output over small window

% for j = 1:length(hbmlp13)
%     dX(j) = hbmlp13(j).intT(2)-hbmlp13(j).intT(1);
%     %find where intT is 10 min
%     min5 = find(hbmlp13(j).intT > 10);
%     for i = 1:length(hbmlp13(j).active)
%         [hbmlp13(j).OutputWindowM(i)] = trapz(hbmlp13(j).intT(1:min5(1)),hbmlp13(j).intM(1:min5(1),hbmlp13(j).active(i)))/(dX(j)*min5(1));
%         [hbmlp13(j).OutputWindowP(i)] = trapz(hbmlp13(j).intT(1:min5(1)),hbmlp13(j).intP(1:min5(1),hbmlp13(j).active(i)))/(dX(j)*min5(1));
%     end
% end
%{
%Look at output from activation to % of peak

for j = 1:length(hbmlp13)
    dX(j) = hbmlp13(j).intT(2)-hbmlp13(j).intT(1);

    for i = 1:length(hbmlp13(j).active)
        begin(i) = find(hbmlp13(j).intT == hbmlp13(j).activenuc(i)); %start at when nucleus becomes active, based off initial activation threshold
        tmp = find(hbmlp13(j).intT > hbmlp13(j).actM(i)); %find when nucleus hits % of peak
        halt(i) = tmp(1);
        timepoints(i) = halt(i)-begin(i);
        if timepoints(i) < 1
            [hbmlp13(j).OutputPercentM(i)] = NaN;
            [hbmlp13(j).OutputPercentP(i)] = NaN;
        else
            [hbmlp13(j).OutputPercentM(i)] = trapz(hbmlp13(j).intT(begin(i):halt(i)),hbmlp13(j).intM(begin(i):halt(i),hbmlp13(j).active(i)))/(dX(j)*(timepoints(i)+1));
            [hbmlp13(j).OutputPercentP(i)] = trapz(hbmlp13(j).intT(begin(i):halt(i)),hbmlp13(j).intP(begin(i):halt(i),hbmlp13(j).active(i)))/(dX(j)*(timepoints(i)+1));
        end
        clear tmp 
    end
end

%Normalize (?) Output

for j = 1:length(hbmlp13)
     maxhb(j) = max(max([hbmlp13(j).intM]));
end
maxhb = max(maxhb);
for j = 1:length(hbmlp13)
    [hbmlp13(j).NormSignalM] = [hbmlp13(j).intM]./maxhb;
end

for j = 1:length(hbmlp13)
    dX(j) = hbmlp13(j).intT(2)-hbmlp13(j).intT(1);
    for i = 1:length(hbmlp13(j).active)
%         [hbmlp13(j).NormOutputM(i)] = trapz(hbmlp13(j).normT,hbmlp13(j).NormSignalM(:,hbmlp13(j).active(i)));%/(dX(j)*100); %integrate using normalized time
        [hbmlp13(j).NormOutputM(i)] = trapz(hbmlp13(j).intT,hbmlp13(j).NormSignalM(:,hbmlp13(j).active(i)))/(dX(j)*100); %integrate using interpolated time and then normalize using total time
    end
end

%}
%% Find distance between max pixels
% flag = [];
% for i = 1:length(hbmlp13)
%     if isempty(hbmlp13(i).active)
%         continue
%     else
%         for j = 1:length(hbmlp13(i).active)
%             for q = 1:length(hbmlp13(i).T)
%                 [hbmlp13(i).dist{j}(q,1)] = pdist([hbmlp13(i).MS2pix(q,:,hbmlp13(i).active(j)); hbmlp13(i).PP7pix(q,:,hbmlp13(i).active(j))]);
%                 if hbmlp13(i).M(q,hbmlp13(i).active(j)) >= threshM & hbmlp13(i).P(q,hbmlp13(i).active(j)) >= threshP
%                     if hbmlp13(i).dist{j}(q,1) < 5
%                         [hbmlp13(i).dist{j}(q,2)] = 1;
%                     else
%                         [hbmlp13(i).dist{j}(q,2)] = 0;
%                         flag(i) = [flag(i) j];
%                     end
%                 else
%                     [hbmlp13(i).dist{j}(q,2)] = NaN;
%                 end
%             end
%         end
% 
%     end
% end


%% Cross-Correlation 
%{
for j = 1:length(hbmlp13)
    for i = 1:length(hbmlp13(j).active)
        [r{j}(:,i),lagspace{j}(:,i)] = xcorr(hbmlp13(j).trueM(:,hbmlp13(j).active(i)),hbmlp13(j).trueP(:,hbmlp13(j).active(i)),'coeff');
        [~,I] = max(r{j}(:,i));
        lag{j}(i) = -1*(lagspace{j}(I,i));
        if lag{j}(i) <= 0
            lag{j}(i) = NaN;
            Erate_corr{j}(i) = NaN;
        else
            Erate_corr{j}(i) = insert/(hbmlp13(j).intT(lag{j}(i)+1)-hbmlp13(j).intT(1));
            if Erate_corr{j}(i) > 10
                Erate_corr{j}(i) = NaN;
            end
        end
        [hbmlp13(j).xcorr(i)] = Erate_corr{j}(i);
    end
end
%}

%% Cross-Correlation Normalized
for j = 1:length(hbmlp13)
    for i = 1:length(hbmlp13(j).active)
        [r{j}(:,i),lagspace{j}(:,i)] = xcorr(hbmlp13(j).normM(:,hbmlp13(j).active(i)),hbmlp13(j).normP(:,hbmlp13(j).active(i)),'coeff');
        [~,I] = max(r{j}(:,i));
        lag{j}(i) = -1*(lagspace{j}(I,i));
        if lag{j}(i) <= 0
            lag{j}(i) = NaN;
            Erate_corr{j}(i) = NaN;
        else
            Erate_corr{j}(i) = insert/(hbmlp13(j).intT(lag{j}(i)+1)-hbmlp13(j).intT(1));
            if Erate_corr{j}(i) > 10
                Erate_corr{j}(i) = NaN;
            end
        end
        [hbmlp13(j).xcorr(i)] = Erate_corr{j}(i);
    end
end

end

