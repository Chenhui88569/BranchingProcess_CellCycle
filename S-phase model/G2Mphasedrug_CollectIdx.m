function indices  = G2Mphasedrug_CollectIdx( AbsolutePhase, obj, str_q)

% Create a cell array to represent your sequence. Replace this with your actual sequence.
states = {'G1_2', 'S', 'G2M', 'G1(2)', 'UR_G2M', 'FR_G2M', 'MS_G2M', 'UR_G1(2)', 'UR_S', 'S(2)', 'G2M(2)'};
  % Number of times the repeating sequence appears

% The repeating sequence
repeatSeq_Treated = {'G1(2)', 'UR_G2M', 'FR_G2M', 'MS_G2M', 'UR_G1(2)', 'UR_S', 'S(2)' , 'G2M(2)'};
repeatSeq_healthy = {'G1' ,'S','G2M'};



% Combine the sequences
fullSeq_Treated = [states, repmat(repeatSeq_Treated, 1, obj.num_RroundsOfTreatedcells-1)];

fullSeq = [fullSeq_Treated,   repmat(repeatSeq_healthy, 1, obj.num_GenOfHealthycells-1 ) ];

if AbsolutePhase== 0
    fullSeq = fullSeq(3:end);
elseif  AbsolutePhase==1
    fullSeq = fullSeq;
elseif  AbsolutePhase== 2
    fullSeq = fullSeq(2:end);
end


% Find the indices containing 'G2M'
indices = find(cellfun(@(x) ~isempty(strfind(x, str_q )), fullSeq ));


end