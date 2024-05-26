% PARSE         parse a vector to sets of points, each of consecutive values.

function [ pairs, idx ] = parse(vec)

if isempty(vec), pairs = []; 
    return
end
if ~isnumeric( vec ) | issparse( vec ) | ~any(size(vec)==1)
    error('input must be a numeric full vector')
end

vec = double( vec(:).' );
[st et] = parsec(vec);
pairs = [st et];
if nargout > 1
    [ ign sidx ] = ismember( st, vec );
    [ ign eidx ] = ismember( et, vec );
    idx = [ sidx eidx ];
end

return

