function [a_G1_gen , a_S_gen, a_G2M_gen, G1_idx_col , S_idx_col,  G2M_idx_col  ] = Cotreatment_CollectIdx_CrossGens(AbsolutePhase, model_obj)
%COTREATMENT_COLLECTIDX Summary of this function goes here
%   Detailed explanation goes here
% Initialize the initial states
state_list = {};

% Fixed states
%state_list = [state_list, 'G1(1)', 'S(1)', 'G2M(1)', 'G1(2)', 'UR_G2M(1)', 'FR_G2M(1)', 'MS_G2M(1)', 'UR_G1(2)', 'UR_S(2)', 'S(2)', 'G2M(2)'];% Repeating states
% Repeating states
for k = 1: model_obj.num_RroundsOfTreatedcells
    currectorder = (k-1)*2+1;
    additional_states = { ...
        ['G1(', num2str(  currectorder  ), ')'], ...
        ['G1_block(', num2str(  currectorder  ), ')'], ...
        ['S(', num2str(  currectorder  ), ')'], ...
        ['UR_S(', num2str(  currectorder ), ')'], ...
        ['FR_S(', num2str(  currectorder  ), ')'], ...
        ['G2M(', num2str(  currectorder  ), ')'], ...
        ['UR_G2M(', num2str(  currectorder  ), ')'] ...
        ['FR_G2M(', num2str(  currectorder  ), ')'] ...
        ['MS_G2M(', num2str(  currectorder  ), ')'] ...
        ['G1(', num2str(  currectorder +1), ')'] ...
        ['G1_block(', num2str(  currectorder  +1), ')'] ...
        ['UR_G1(', num2str(  currectorder  +1), ')'] ...
        ['S(', num2str(  currectorder  +1), ')'] ...
        ['UR_S(', num2str(  currectorder  +1), ')'] ...
        ['FR_S(', num2str(  currectorder  +1), ')'] ...
        ['G2M(', num2str(  currectorder  +1), ')'] ...
        ['UR_G2M(', num2str(  currectorder  +1), ')'] ...
        ['FR_G2M(', num2str(  currectorder  +1), ')'] ...
        ['MS_G2M(', num2str(  currectorder  +1), ')'] ...
        };
    
    state_list = [state_list, additional_states];
end

if model_obj.num_GenOfHealthycells-1>0
    for k = 1:model_obj.num_GenOfHealthycells-1
        order = currectorder  +2 +k;
        additional_states = { ...
            ['G1(', num2str( order ), ')'], ...
            ['S(', num2str( order  ), ')'] ...
            ['G2M(', num2str( order  ), ')'] ...
            };
        
        state_list = [state_list, additional_states];
    end
else
    order =  currectorder  + 1;
end

if     AbsolutePhase == 0
    state_list  =  state_list(1:end);
elseif     AbsolutePhase == 1
    state_list  =  state_list(2:end);
elseif     AbsolutePhase == 2
    state_list  =  state_list(4:end);
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

