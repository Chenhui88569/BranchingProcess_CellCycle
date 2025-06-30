function [a_G1_gen , a_S_gen, a_G2M_gen, G1_idx_col , S_idx_col,  G2M_idx_col  ] = SequentialTreatment_CollectIdx_CrossGens(AbsolutePhase, obj)
%COTREATMENT_COLLECTIDX Summary of this function goes here
%   Detailed explanation goes here
% Initialize the initial states
filename = 'transition_matrix.xlsx';
if obj.TreatmentEffectOnPhases_1thDrug == 3 && obj.TreatmentEffectOnPhases_2thDrug == 2 % G2M phase first
    
    % Read the entire table (assuming 'data.xlsx' is your Excel file name)
    T = readtable(filename, 'Sheet' ,'G2M phase first' ,'VariableNamingRule','preserve');
    % Extract the headers (variable names)
    headers = T.Properties.VariableNames;
    state_list = headers;
    if     AbsolutePhase == 0
        state_list  =  state_list(2:end);
    elseif     AbsolutePhase == 1
        state_list  =  state_list(1:end);
    elseif     AbsolutePhase == 2
        state_list  =  state_list(2:end);
    end
    
    
elseif  obj.TreatmentEffectOnPhases_1thDrug == 2 && obj.TreatmentEffectOnPhases_2thDrug == 3 % S first
    T = readtable(filename, 'Sheet' ,'S phase first' ,'VariableNamingRule','preserve');
    % Extract the headers (variable names)
    
    headers = T.Properties.VariableNames;
    headers  = headers(2:end);
    % Function to increment numbers within the parentheses
   % incrementFunc = @(match) ['(' num2str(str2double(match{1}(2:end-1)) + 1) ')'];
    
    % Process each header entry
    newHeader = cell(size(  headers  ));
    for i = 1:numel(  headers  )
        newHeader{i} = regexprep(headers{i}, '\((\d+)\)', '\(${num2str(str2double($1) + 1)}\)');
    end
    state_list = newHeader;
    if     AbsolutePhase == 0
        state_list  =   [ 'G2M_2(1)', state_list];
        state_list  =  state_list(4:end);
    elseif     AbsolutePhase == 1
        state_list  =  state_list(1:end);
    elseif     AbsolutePhase == 2
        state_list  =  state_list(3:end);
    end


end



a_G1_gen = cell(  order, 1 );
a_S_gen =  cell(  order, 1 );
a_G2M_gen =  cell(  order, 1 );

% Loop through the state_list to filter the names and extract numbers
for i = 1:length(state_list)
    state = state_list{i};
    if contains(state, 'G1')
        num = regexp(state, '\((\d+)\)', 'tokens');
        num = str2double(num{1}{1});
        a_G1_gen{num} = [a_G1_gen{num}, i];
    elseif contains(state, 'S')
        num = regexp(state, '\((\d+)\)', 'tokens');
        num = str2double(num{1}{1});
        a_S_gen{num} = [a_S_gen{num}, i];
    elseif contains(state, 'G2M')
        num = regexp(state, '\((\d+)\)', 'tokens');
        num = str2double(num{1}{1});
        a_G2M_gen{num} = [a_G2M_gen{num}, i];
    end
end

G1_idx_col  = find(cellfun(@(x) ~isempty(strfind(x, 'G1' )), state_list  )); 
S_idx_col= find(cellfun(@(x) ~isempty(strfind(x, 'S' )), state_list  )); 
G2M_idx_col= find(cellfun(@(x) ~isempty(strfind(x, 'G2M' )), state_list  )); 
end

