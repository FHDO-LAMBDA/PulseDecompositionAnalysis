function [singleBeats_processed] = createSingleBeats(PPG,samplingFrequency,beatIndices,beatIndicesEnsembleBeat)
%
%
% add inputs and outputs here!
%
%

%% Preprocessing
% turn ppg_signal into 1xn array
if(size(PPG,1)>size(PPG,2))
    PPG = PPG';
end

% filter as done with the ensemble data
cutOff_lower = 0.2; % lower cutoff frequency
cutOff_upper = 8; % upper cutoff frequency
order = 5; % filter order
[coeff_b,coeff_a] = butter(order, [cutOff_lower cutOff_upper]/(samplingFrequency/2), 'bandpass'); %create bandpass to filter the pixel traces
PPG = filtfilt(coeff_b, coeff_a, PPG);

%% Extract single beats
segmentLength = diff(beatIndices); % calculate beat to beat differences
segmentLength = ceil(median(segmentLength)); % get median segment length
beatInterval_before = round(0.45*segmentLength); % beat interval before detection point
beatInterval_after = segmentLength; % beat interval after detection point
singleBeats = arrayfun(@(x) PPG(x-beatInterval_before:x+beatInterval_after),beatIndicesEnsembleBeat,'UniformOutput', false); % cut PPG into single beats
singleBeats_processed = cell(size(singleBeats)); % initialize cell array for processed beats

for beatNumber = 1:size(singleBeats,1)
    %% Find start and end of each beat
    currentBeat = singleBeats{beatNumber};
    
    % wenn ich hier keine start und end-Punkte finde, muss ich ja eigentlich
    % den ganzen beat raussschmeißen. Das will ich aber eigentlich nicht, weil
    % jeder dieser beats in mein ensemble eingegangen ist und ich jetzt aus
    % allen einzel beats jeweils die Parameter berechne und dann z.B. über alle
    % Einzelparameter mitteln will statt alle beats zu mitteln
    
    % ich brauche die Filterung aber, um die Start und Endpunkte so nutzen
    % zu können wie sebastian! ich kann nur nach der minima detektion
    % entscheiden, ob ich ungefilterte daten nehme oder gefilterte

    % get initial minimum (start)
    startingSegment = currentBeat(1:beatInterval_before); % get frst part of the beat
    [minima,minimaIndices] = findpeaks(-startingSegment); % get minima
    if(~isempty(minima))
        index = numel(minimaIndices);%take last minimum before slope
        beatStartIndex = minimaIndices(index);
    else
        beatStartIndex = 1; % einfach frühstmöglichen Punkt nehmen --> Sebastian hat ja diese schläge auch alle genutzt
    end
    clear minima minimaIndices % needed as the number of detected minima can vary between different beats
    
    % get ending minimum (end)
    endingSegment = currentBeat(beatInterval_before:end); % get second part of the beat
    [minima,minimaIndices] = findpeaks(-endingSegment); % get minima
    if(~isempty(minima))
        [~,index] = min(-minima);%take lowest minimum
        beatStopIndex = beatInterval_before + minimaIndices(index);
    else
        beatStopIndex = length(currentBeat); % einfach letzten Punkt nehmen
    end
    clear minima minimaIndices % needed as the number of detected minima can vary between different beats    
    
    %% Detrend beat
    % interpolate a straight line between the detected minima and
    % extrapolate that line to the length of the total segment, then
    % subtract that from the segment
    trenddata = interp1([beatStartIndex beatStopIndex], [currentBeat(beatStartIndex) currentBeat(beatStopIndex)],1:length(currentBeat),'linear','extrap');%calculate straight line from first to last point
    currentBeat = currentBeat - trenddata;%do detrending by removing the straight line
    
    %% Shorten beat
    currentBeat = currentBeat(beatStartIndex:beatStopIndex);
    singleBeats_processed{beatNumber} = currentBeat;
    clear currentBeat startingSegment endingSegment trenddata
end

end