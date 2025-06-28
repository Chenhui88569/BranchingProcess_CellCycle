function parts = splitYamlString(str)
% Extract LaTeX-style parameters from a quoted string list like:
% '"$q_1$","$EC_{50,d}$","$n_{G2M}$"'

    % Use regular expression to extract text within double quotes
    parts = regexp(str, '"([^"]+)"', 'tokens');
    
    % Flatten the cell array of cells
    parts = [parts{:}];
end
