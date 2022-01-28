% get diff sets among over than two vector sets
% varargin{1}: shared elements among all sets
function diffsets_tot = getdiffsets(varargin)
ret = varargin{1};
diffsets_tot = [];
for k = 2:nargin
    [diffsets,~] = setdiff(varargin{k},ret,'stable');
    diffsets_tot = [diffsets_tot;diffsets((1:100))];
end

end