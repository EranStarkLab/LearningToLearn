% LtL_MAKE_EQUAL_BINS       assign indices to data s.t. equal sized bins.
%
% call                  [ SVAL EVAL IDX N ] = LtL_MAKE_EQUAL_BINS( X, N )
% 
% gets                  X           data
%                       N           number of bins
%                       DIF         {1}, discrete interger flag
%
% return                SVAL        lowest value in each bin (length of N)
%                       EVAL        highest value in each bin (")
%                       IDX         bin numbers (same size as X)
%                       N           number of bins de-facto

% 08-nov-04 ES


function [ sval, eval, idx, n ] = LtL_make_equal_bins( x, n, dif )

nargs = nargin;
if nargs < 2 || isempty( x ) || isempty( n ), error( '2 arguments' ), end
if nargs < 3 || isempty( dif )
    if all( x( : ) == round( x( : ) ) )
        dif = 1;
    else
        dif = 0;
    end
end

% prepare
sx = size( x );
if any( sx == 1 )                               % single vector input
    x = x( : );
    sx = size( x );
end
nnans = ~isnan( x );
if sum( sum( ~nnans ) ) && dif == 1
    if sx( 2 ) == 1
        x( ~nnans ) = []; % remove all NaNs
    else
        error( 'not supported' )
    end
end
lx = sum( nnans );

% determine edges
if dif == 1
    % initialize
    N = n;
    sval = zeros( n( sx( 2 ) ), sx( 2 ) );
    eval = sval;
    for c = 1 : sx( 2 )
        %eii = [];
        [ h v ] = uhist( x( :, c ) );                        % make a histogram of discrete values
        nh = length( h );                                   % number of different values
        f = cumsum( h / lx( c ) );                          % cumulative frequency
        n( c ) = min( N, ceil( lx( c ) / max( h ) ) );      % upper bound for number of bins
        eii = zeros( n( c ), 1 );
        for i = 1 : ( n( c ) - 1 )
            bin_end = find( f >= i / n( c ) , 1 );
            if isempty( bin_end ) || bin_end == nh
                bin_end = find( f < i / n( c ) , 1, 'last' );
            end
            eii( i ) = bin_end;
        end
        eii = eii( : );
        eii( n( c ) ) = nh;
        sii = [ 1; eii( 1 : n( c ) - 1 ) + 1 ];
        sval( 1 : n( c ), c ) = v( sii )';
        eval( 1 : n( c ), c ) = v( eii )';
    end
else
    [ xs ix ] = sort( x );
    bsize = lx / n;
    si = cumsum( [ ones( 1, sx( 2 ) ); ones( n - 1, 1 ) * bsize ], 1 );
    ei = cumsum( ones( n, 1 ) * bsize, 1 );
    ios = ones( n, 1 ) * [ 0 : sx( 1 ) : ( sx( 1 ) * ( sx( 2 ) - 1 ) ) ];
    si = floor( si + ios );
    ei = ceil( ei + ios );
    si( si < 1 ) = 1;
    ei( ei > lx ) = lx;
    sval = xs( si );
    eval = xs( ei );
    
    % 20-jun-14 hack for identicals 
    rmv = find( diff( sval ) == 0 );
    if ~isempty( rmv )
        sval( rmv ) = [];
        eval( rmv ) = [];
    end
    rmv = find( diff( eval ) == 0 );
    if ~isempty( rmv )
        sval( rmv ) = [];
        eval( rmv ) = [];
    end
end

% post-proc
if dif == 1
    sevals = sortranges( [ sval eval ], 0 );
    sevals  = unique( sevals, 'rows' );
    sval = sevals( :, 1 );
    eval = sevals( :, 2 );
end

% actually bin
idx = zeros( sx );
if nargout > 2
    if dif
        for c = 1 : sx( 2 )
            for i = 1 : n( c )
                idx( ( x( :, c ) >= sval( i, c ) ) & ( x( :, c ) <= eval( i, c ) ), c ) = i;
            end
        end
    else
        for i = 1 : n
            for c = 1 : sx( 2 )
                idx( ix( si( i, c ) : ei( i, c ) ), c ) = i;
            end
        end
    end
end

return

