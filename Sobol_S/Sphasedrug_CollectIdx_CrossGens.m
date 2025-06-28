function [a_G1, a_S, a_G2M]  = Sphasedrug_CollectIdx_CrossGens( AbsolutePhase, model_obj)

% Initialize the initial states
state_list = {};

% Fixed states
%state_list = [state_list, 'G1(1)', 'S(1)', 'G2M(1)', 'G1(2)', 'UR_G2M(1)', 'FR_G2M(1)', 'MS_G2M(1)', 'UR_G1(2)', 'UR_S(2)', 'S(2)', 'G2M(2)'];% Repeating states
% Repeating states
state_list  =   [ 'G1M_2(1)', state_list];
for k = 2: model_obj.num_RroundsOfTreatedcells+1
    additional_states = { ...
        ['G1(', num2str(k), ')'], ...
        ['G1_block(', num2str(k), ')'], ...
        ['S(', num2str(k), ')'], ...
        ['G2M(', num2str(k), ')'], ...
        ['UR_S(', num2str(k), ')'], ...
        ['FR_S(', num2str(k), ')'], ...
        ['UR_G2M(', num2str(k), ')'] ...
        };
    
    state_list = [state_list, additional_states];
end

if model_obj.num_GenOfHealthycells-1>0
    for k = 1:model_obj.num_GenOfHealthycells-1
        order = model_obj.num_RroundsOfTreatedcells+k+1;
        additional_states = { ...
            ['G1(', num2str( order ), ')'], ...
            ['S(', num2str( order  ), ')'] ...
            ['G2M(', num2str( order  ), ')'] ...
            };
        
        state_list = [state_list, additional_states];
    end
else
    order = model_obj.num_RroundsOfTreatedcells+1;
end

if     AbsolutePhase == 0
    state_list  =  state_list(1:end);
elseif     AbsolutePhase == 1
    state_list  =  state_list(2:end);
elseif     AbsolutePhase == 2
    state_list  =  state_list(4:end);
end

a_G1 = cell(  order, 1 );
a_S =  cell(  order, 1 );
a_G2M =  cell(  order, 1 );

% Loop through the state_list to filter the names and extract numbers
for i = 1:length(state_list)
    state = state_list{i};
    if contains(state, 'G1')
        num = regexp(state, '\((\d+)\)', 'tokens');
        num = str2double(num{1}{1});
        a_G1{num} = [a_G1{num}, i];
    elseif contains(state, 'S')
        num = regexp(state, '\((\d+)\)', 'tokens');
        num = str2double(num{1}{1});
        a_S{num} = [a_S{num}, i];
    elseif contains(state, 'G2M')
        num = regexp(state, '\((\d+)\)', 'tokens');
        num = str2double(num{1}{1});
        a_G2M{num} = [a_G2M{num}, i];
    end
end

% % Convert cell arrays to numeric arrays
% a_G1 = cell2mat(a_G1);
% a_S = cell2mat(a_S);
% a_G2M = cell2mat(a_G2M);
%
% % Display the matrices
% disp('Numbers for G1:');
% disp(a_G1);
% disp('Numbers for S:');
% disp(a_S);
% disp('Numbers for G2M:');
% disp(a_G2M);
%
%
% disp(state_list)



end