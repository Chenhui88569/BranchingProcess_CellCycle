function indices  = Sphasedrug_CollectIdx( AbsolutePhase, obj, str_q)

% Create a cell array to represent your sequence. Replace this with your actual sequence.
states =  {'G1_2', 'G1_block', 'S', 'G2M', 'UR_S', 'FR_S', 'UR_G2M'};
  % Number of times the repeating sequence appears

% The repeating sequence
repeatSeq_Treated = {'G1_2', 'G1_block', 'S', 'G2M', 'UR_S', 'FR_S', 'UR_G2M'};
repeatSeq_healthy = {'G1' ,'S','G2M'};



% Combine the sequences
fullSeq_Treated = [states, repmat(repeatSeq_Treated, 1, obj.num_RroundsOfTreatedcells-1)];

fullSeq = [fullSeq_Treated,   repmat(repeatSeq_healthy, 1, obj.num_GenOfHealthycells-1 ) ];

if AbsolutePhase== 0
    fullSeq =  ['G2M', states];
elseif  AbsolutePhase==1
    fullSeq = fullSeq;
elseif  AbsolutePhase== 2
    fullSeq = fullSeq(3:end);
end


% Find the indices containing 'G2M'
indices = find(cellfun(@(x) ~isempty(strfind(x, str_q )), fullSeq ));


%disp('Indices containing "G2M":');
%disp(indices);
end