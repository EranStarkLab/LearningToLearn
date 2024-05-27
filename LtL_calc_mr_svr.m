% LtL_calc_mr_svr               SVR/QR (regression) or SVM (classification)
%
% does                          cross-validated regresssion/classification
%
% call                          [ R2, res, fig ] = LtL_calc_mr_svr( x, y )
%
% returns                       R2 (svr/qr) or PC (svm)
%
% calls                         LtL_ParseArgPairs
%                               svmtrain, svmpredict
%                               LtL_make_equal_bins
%                               LtL_mutinf, LtL_calc_bias_PT

% 18-dec-22 ES

% Last update:
% 26-may-24 AL

% procedure for SVM/SVR
% (1) divide data points (rows) randomly into folds, e.g., 10 folds
% (2) for every fold, train the model on the e.g., 90%, then predict on the 10%
% (3) in every fold, first run a grid search, then use those parameters to
% train the final model
% (4) after repeating for all folds, quantify using accuracy/AUC/MI
% (classification) or R2 (regression)
% (5) to compute the contribution of every feature (columns), remove that
% feature and redo everything

function [ R2, res, fig ] = LtL_calc_mr_svr( x, y, varargin )

% defaults
preprocess_DFLT                 = 'scale';
model_DFLT                      = 'svr';
svm_grid_search_DFLT            = 1;
nt_grid_search_DFLT             = NaN;
svm_kernel_DFLT                 = 2; %0/1/{2}/3 (linear/polynomial/RBF/sigmoid)
nfold_DFLT                      = 2;
nbins_DFLT                      = 4;

% initialize output
R2                              = [];
res                             = [];
fig                             = [];

% arguments
nargs                           = nargin;
if nargs < 2 || isempty( x ) || isempty( y )
    return
end
[ nt, nc ]                      = size( y );
nt2                             = size( x, 1 );
if nt ~= nt2
    error( 'input size mismatch: must have same number of samples (rows) in x and y' )
end
if nc ~= 1
    error( 'y can have only 1 column' )
end

[ preprocess, model, svm_grid_search, nt_grid_search, svm_kernel ...
    , nfold, nbins ...
    , graphics, verbose ]       = LtL_ParseArgPairs(...
    { 'preprocess', 'model', 'svm_grid_search', 'nt_grid_search', 'svm_kernel' ...
    , 'nfold', 'nbins' ...
    , 'graphics', 'verbose' }...
    , { preprocess_DFLT, model_DFLT, svm_grid_search_DFLT, nt_grid_search_DFLT, svm_kernel_DFLT ...
    , nfold_DFLT, nbins_DFLT ...
    , 0, 0 }...
    , varargin{ : } );

if ~ismember( preprocess, { 'scale', 'zs', 'standardize', 'rank', 'none', 'null' } )
    error( 'unsupported preprocess' )
end
if ~ismember( model, { 'svr', 'qr', 'svm' } )
    error( 'unsupported model' )
end
if ~ismember( svm_grid_search, [ 0 1 ] )
    error( 'unsupported svm_grid_search' )
end
if isnan( nt_grid_search )
    nt_grid_search           	= nt;
end

if ~ismember( svm_kernel, [ 0 1 2 3 ] )
    error( 'unsupported svm_kernel' )
end
if isnumeric( nfold ) && length( nfold ) ~= 1 && ~isequal( nfold, 'all' )
    error( 'unsupported nfold' )
end
if isa( nfold, 'char' ) && isequal( nfold, 'all' )
    nfold                       = nt;
    testTrial                   = [];
elseif round( nfold ) ~= nfold
    testTrial                   = nfold * 1000;
    nfold                       = 1;
end
if isnan( nbins )
    nbins                       = 2 ^ ( nextpow2( sqrt( nt / 10 ) ) - 1 );
end

% preprocess the data
switch preprocess
    case { 'scale' }
        x                       = scale_local( x );
    case { 'zs', 'standardize' }
        x                       = zs_local( x );
    case { 'rank' }
        x                       = rankcols_local( x );
        y                       = rankcols_local( y );
end


% prepare for SVM
if strcmp( model, 'svr' )
    svm_str                     = sprintf( '-d %d, -s %d, -t %d', 3, 3, svm_kernel );
    if svm_grid_search                          % grid search each parameter
        nfold_grid_search       = 2;
        ce                      = 2 .^( -5 : 2 : 15 );
        ge                      = 2 .^( -15 : 2 : 3 );
        nc                      = length( ce );
        ng                      = length( ge );
        [ ~, tmp ]              = sort( rand( nt, 1 ) );
        ts                      = sort( tmp( 1 : min( nt, nt_grid_search ) ) );  % select nt_grid_search random trials
        Xgrid                   = x( ts, : );
        Ygrid                   = y( ts, : );
        if verbose
            fprintf( 1, 'Grid searching parameter...' )
        end
        f                       = NaN( nc, ng );
        for ci                  = 1 : nc
            if verbose
                fprintf( 1, '*' )
            end
            for gi              = 1 : ng
                str             = sprintf( '-v %d, -c %0.3g, -g %0.3g, %s'...
                    , nfold_grid_search, ce( ci ), ge( gi ), svm_str );
                f( ci, gi )     = svmtrain( Ygrid, Xgrid, str );
            end
        end
        if verbose
            fprintf( 1, '\n' )
        end
        [ i, j ]              	= find( f == min( min( f ) ) );
        c                       = ce( i( 1 ) );
        g                       = ge( j( 1 ) );
        %figure
        %imagesc( log2( ge ), log2( ce ), f( :, :, m ) ), axis xy
        %xlabel( 'log( gamma )' ), ylabel( 'log( cost )' )
        %set( gca, 'xtick', log2( ge ), 'ytick', log2( ce ) )
        %title( sprintf( 'm = %d', m ) )
    else
        c                       = 5;%1;
        g                       = 2 ^ -7;
    end
    res.c                       = c;
    res.g                       = g;
    svm_params                  = sprintf( '-c %0.3g, -g %0.3g, %s'...
        , c, g, svm_str );
    % 3, 3 are svmtrain parameters: 3rd degree polynomial; epsilon-SVR
    % we don't learn the coef0 of a polynomial kernel
elseif strcmp( model, 'svm' )
    svm_scale                   = 0;
    [ c, g ]                    = svmgrid( x, y, svm_scale, svm_kernel, verbose );
    svm_params                  = sprintf( '-c %0.3g, -g %0.3g, -t %d'...
        , c, g, svm_kernel );
end

% divide trials into nfold disjoint sets
if nfold ~= 1 && nfold ~= nt
    [ ~, bb ]                   = sort( rand( nt, 1 ) );
    div                         = [ 0 ceil( nt / nfold : nt / nfold : nt ) ];
    tidxmat                     = true( nt, nfold );
    for ni                      = 1 : nfold
        idx                     = bb( ( div( ni ) + 1 ) : div( ni + 1 ) );
        tidxmat( idx, ni )      = 0;
    end
end


% regress/classify using n-fold cross-validation
if verbose
    switch model
        case { 'qr', 'svr' }
            fprintf( 1, 'Regressing (%d-fold cross-validation)', nfold )
        case 'svm'
            fprintf( 1, 'Classifying (%d-fold cross-validation)', nfold )
    end
end
YtestAccum                      = NaN( nt, 1 );
YhatAccum                       = NaN( nt, 1 );
for ni = 1 : nfold
    if verbose
        fprintf( 1, '.' )
    end
    % divide trials into train and test sets
    if nfold == 1                       % train and test are the same
        trainidx                = true( nt, 1 );
        testidx                 = trainidx;
        if exist( 'testTrial', 'var' ) && testTrial <= nt
            trainidx( testTrial ) = 0;  % test a specific trial
            testidx             = ~trainidx;
        end
    elseif nfold == nt                  % leave-one-out sequentially
        trainidx                = true( nt, 1 );
        trainidx( ni )          = 0;
        testidx                 = ~trainidx;
    elseif nfold > 0                    % nfold disjoint groups
        trainidx                = tidxmat( :, ni );
        testidx                 = ~trainidx;
    else                                % select subset randomly at each iteration
        trainidx = ~true( NT, 1 );
        [ ~, bb ]               = sort( rand( nt, 1 ) );
        oidx                    = bb( 1 : round( nt * ( 1 - 1 / -nfold ) ) );
        trainidx( oidx )    	= 1;
        testidx                 = ~tidx;
    end
    % training set
    Xtrain                      = x( trainidx, : );
    Ytrain                      = y( trainidx, : );
    if strcmp( preprocess, 'standardize' )
        Xtrain                  = zs_local( Xtrain );
        Ytrain                  = zs_local( Ytrain );
    end
    % there are two differences between the 3 models:
    % (1) MISO/MIMO - all are MISO, but the linear solves simultaneuosly
    % (2) test trials - linear and SVM are trial-by-trial, SCG is concatenated
    Xtest                       = x( testidx, : );
    Ytest                       = y( testidx, : );
    
    % exclude sparse coding
    ridx                        = sum( Xtrain.^2 ) == 0;
    if sum( ridx ) > 0
        fprintf( '%s: %d''th fold: %d/%d columns empty!!\n'...
            , upper( mfilename ), ni, sum( ridx ), length( ridx ) )
        Xtrain( :, ridx )     	= [];
        Xtest( :, ridx )        = [];
    end
    
    % regression and reconstruction
    ntest                       = sum( testidx );
    Yhat                        = NaN( ntest, 1 );
    switch model
        case { 'svr', 'svm' } % support vector machine
            % since the current implementation of support vector regression
            % does not support a MIMO system, we regress M separate MISO
            % systems (for the linear system it is the same)
            % train a model for the m'th parameter
            svm_model           = svmtrain( Ytrain, Xtrain, svm_params );
            % and test it
            Yhat                = NaN( ntest, 1 );
            for t               = 1 : ntest
                Yhat( t )       = svmpredict( Ytest( t ), Xtest( t, : ), svm_model );
            end
        case 'qr'   % wiener filter
            % the linear model is		Y = X*A
            % the solution is		    A = inv(X'*X)*X'*Y
            % the reconstruction is		Yhat = X*A;
            % solve using the numerically more stable QR decomposition:
            [ Q, R ]            = qr( Xtrain, 0 );
            A                   = R \ ( Q' * Ytrain );
            % reconstruct each test trial separately
            for t               = 1 : ntest
                Yhat( t )       = Xtest( t, : ) * A;
            end
    end
    % regression statistics
    YtestAccum( testidx )       = Ytest;
    YhatAccum( testidx )        = Yhat;
end

switch model
    case { 'qr', 'svr' }
        [ R2, mi, rmat ]        = regstats( YhatAccum, YtestAccum, nbins );
    case 'svm'
        [ R2, mi, rmat ]        = classstats( YhatAccum, YtestAccum );
end
res.model                       = model;
res.R2                          = R2;                                       % PC for a classification model
res.mi                          = mi;
res.rmat                        = rmat;
res.ytest                       = YtestAccum;
res.yhat                        = YhatAccum;

if ~graphics
    return
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% classstats          classification statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ pc, mi, cmat ] = classstats( yhat, y )

% PC statistics
pc                              = sum( yhat == y ) / length( yhat );

% classification matrix
cls                             = unique( [ yhat; y ] );
if cls( 1 ) == 0
    cls                         = cls + 1;
    yhat                        = yhat + 1;
    y                           = y + 1;
end
C                               = length( cls );
N                               = length( yhat );
cmat                            = zeros( C, C ); % rows - correct, cols - reconstructed
for i                           = 1 : N
    if ismember( yhat( i ), cls )
        cmat( y( i ), yhat( i ) ) = cmat( y( i ), yhat( i ) ) + 1;
    end
end

% mutual information
mi                              = LtL_mutinf( cmat ) - LtL_calc_bias_PT( cmat, 2 );

return % classstats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% regstats          regression statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ r2, mi, rmat ] = regstats( yhat, y, C )

% R2 statistics
r2                              = calc_pearson_local( yhat( : ), y( : ) ) .^ 2;
%r2                              = calc_spearman( yhat( : ), y( : ) ) .^ 2;

% reconstruction matrix
N                               = size( y, 1 );
[ ~, ~, clas0 ]                 = LtL_make_equal_bins( y, C );
[ ~, ~, clas ]                  = LtL_make_equal_bins( yhat, C );
rmat                            = zeros( C, C ); % rows - correct, cols - reconstructed
for i                           = 1 : N
    rmat( clas0( i ), clas( i ) ) = rmat( clas0( i ), clas( i ) ) + 1;
end

% mutual information
mi                              = LtL_mutinf( rmat ) - LtL_calc_bias_PT( rmat, 2 );

return % regstats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% svmgrid           find data-dependent C and gamma for SVM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% call              [ c, g, f ] = svmgrid( x, y, scale );
% gets              x, y        activity matrix (NxP) and tags vector(Nx1)
%                   scale       {0}; 1 means to scale columns of x to the range [0-1]
% returns           c, g        cost and gamma (RBF SD) parameters of SVM
%                   f           the 2-fold accuracy per {c,g} set tested
%
% NOTE              this search is far from perfect - it could be improved
%                   by recursion
function [ c, g, f ] = svmgrid( x, y, scale, kernel, verbose )

if nargin < 3 || isempty( scale )
    scale                       = 0;
end
if scale
    num                         = x - ones( size( x, 1 ), 1 ) * min( x );
    den                         = ones( size( x, 1 ), 1 ) * max( x ) - ones( size( x, 1 ), 1 ) * min( x );
    x                           =  num ./ den;
end

nfold                           = 2;
if kernel == 0
    ge                          = 1;             % linear kernel
    ce                          = -10 : 1 : 3;
else
    ge                          = -15 : 2 : 3;   % polynomial, RBF, and sigmoid kernels
    ce                          = -5 : 2 : 15;
end

f                               = NaN( length( ce ), length( ge ) );
for ci                          = 1 : length( ce )
    if verbose
        fprintf( 1, '#' )
    end
    for gi                      = 1 : length( ge )
        c                       = 2 ^ ce( ci );
        g                       = 2 ^ ge( gi );
        if kernel == 0
            str                 = sprintf( '-v %d, -c %0.3g, -t %d', nfold, c, kernel );
        else
            str                 = sprintf( '-v %d, -c %0.3g, -g %0.3g, -t %d', nfold, c, g, kernel );
        end
        f( ci, gi )             = svmtrain( y, x, str );
    end
end
[ i, j ]                        = find( f == max( max( f ) ) );
c                               = 2 ^ ce( i( 1 ) );
g                               = 2 ^ ge( j( 1 ) );
if verbose
    fprintf( 1, '\n' )
end

return  %LtL_calc_mr_svr

%------------------------------------------------------------------------
% y = scale_local( x )
% scale to the range 0-1
%------------------------------------------------------------------------%

function y = scale_local(x)

d = ones( size( x, 1 ), 1 );
num = x - d * min( x );
den = d * max( x ) - d * min( x );
y =  num ./ den;
y( :, sum( den ) == 0 ) = 0;

return  %scale_local
%------------------------------------------------------------------------%


%------------------------------------------------------------------------%
% [ y, T ] = rankcols_local( x )
% rank matrix column-wise.
%------------------------------------------------------------------------%
% x       data matrix
function [ y, T ] = rankcols_local( x )

[ m, n ]            = size( x );
if m == 1
    x               = x';
    m               = n;
    n               = 1;
end
nans                = isnan( x );
ridx                = m : -1 : 1;
cidx                = 1 : n;
[ x, idx ]          = sort( [ x x( ridx, : ) ] );
[ ~, ranks ]        = sort( idx );
y                   = ( ranks( :, cidx ) + ranks( ridx, cidx + n ) ) / 2;
y( nans )           = NaN;

if nargout > 1
    T               = zeros( 1, n );
    for i           = 1 : n
        t           = diff( parse( find( diff( x( :, i ) ) == 0 ) ), [], 2 ) + 2; 
        T( i )      = sum( ( t - 1 ) .* t .* ( t + 1 ) );
    end
end

return   %rankcols_local
%------------------------------------------------------------------------%


%------------------------------------------------------------------------%
% [ z, m, s, n ] = zs_local( x, tvec, flag )
% z-score for matrix columns.
%------------------------------------------------------------------------%
% x       observations in rows, variables in columns
% TVEC    optional column vector, same for all columns in
%         X. if supplied, each column is standardized in
%         groups, according to the pointer vector TVEC.
% FLAG    {0}; 1 normalizes by N and thus allows
%         computation of correlations
%
% calls             ZSC (mex).

function [ z, m, s, n ] = zs_local( x, tvec, flag )

nargs                       = nargin;
nout                        = nargout;
if nargs < 1 || isempty( x ) || ndims( x ) ~= 2
    error( 'matrix data' )
end
if nargs < 2 || isempty( tvec )
    tvec                    = [];
end
if nargs < 3
    flag                    = 0;
end
[ n, ncols ]                = size( x );
z                           = zeros( n, ncols );
nans                        = isnan( x );

if isempty( tvec )
    
    % standardize each column as a whole
    nanf                    = sum( sum( nans ) );
    if nout <= 1 && ~nanf
        % call a fast routine
        z                   = zsc( double( x ), flag );
    else
        % otherwise, do everthing local
        d                   = ones( n, 1 );
        if nanf
            s               = nanstd( x );
            m               = nanmean( x );
            y               = x - d * m;
        else
            m               = sum( x, 1 ) / n;
            y               = x - d * m;
            s               = sqrt( sum( y .* y, 1 ) / ( n - ~flag ) );
        end
        clear x
        zidx                = s == 0;
        den                 = d * s;
        z( :, ~zidx )       = y( :, ~zidx ) ./ den( :, ~zidx );
    end
    
else
    
    % standardize each column by parts, according to pointer vector tvec
    tvec                    = tvec( : ).';
    if n ~= length( tvec )
        error( 'input size mismatch' )
    end
    uidx                    = unique( tvec );
    
    % go over subsets
    for i                   = uidx( : ).'
        
        idx                 = tvec == i;
        n( i, : )           = sum( idx );
        nanf                = sum( sum( nans( idx ) ) );
        if nout <= 1 && ~nanf
            z( idx, : )     = zsc( x( idx, : ), flag );
        else
            di              = ones( n( i ), 1 );
            xi              = x( idx, : );
            if nanf
                si          = nanstd( xi );
                mi          = nanmean( xi );
                yi          = xi - di * mi;
            else
                mi          = sum( xi, 1 ) / n( i );
                yi          = xi - di * mi;
                si          = sqrt( sum( yi .* yi, 1 ) / ( n( i ) - ~flag ) );
            end
            clear xi
            zidx            = si == 0;
            den             = di * si;
            z( idx, ~zidx ) = yi( :, ~zidx ) ./ den( :, ~zidx );
            m( i, : )       = mi;
            s( i, : )       = si;
        end
        
    end
    
end

return  %zs_local
%------------------------------------------------------------------------%


%------------------------------------------------------------------------
% cc = calc_pearson_local( x, y )
% compute Pearson's correlation coeffiecient
%------------------------------------------------------------------------%
% x       data matrix (2 columns)
% y       optional; then computes the CC between columns of x and y
function cc = calc_pearson_local( x, y )

if isempty( x ) && isempty( y )
    cc                      = [];
    return
end
[ m, n ]                    = size( x );
if m < 2
    cc                      = NaN * ones( 1, n );
    return
end
d                           = ones( m, 1 );

if nargin == 1 || isempty( y )           % correlation matrix (all possible x's column-pairs)
    
    v                       = ~isnan( x );
    vm                      = sum( v );
    y                       = x - d * ( nansum( x .* v ) ./ vm );
    s                       = sqrt( nansum( y .* y .* v ) );
    z                       = y ./ ( d * s );
    z( y == 0 )             = 0;
    z( isnan( z ) )         = 0;
    cc                      = z' * z;
    if n == 2                           % correlation between columns (scalar)
        cc                  = cc( 1, 2 );
    end
    
elseif ~isequal( [ m n ], size( y ) )
    error( 'input size mismatch' )
    
else                                    % correlation vector (column pairs of x,y)
    v                   = ~isnan( x ) & ~isnan( y ); 
    vm                  = sum( v );
    x                   = x - d * ( nansum( x .* v ) ./ vm );
    y                   = y - d * ( nansum( y .* v ) ./ vm );
    num                 = nansum( x .* y .* v );
    den                 = sqrt( nansum( x .* x .* v ) .* nansum( y .* y .* v) );
    zidx                = den == 0;
    cc                  = zeros( 1, n );
    cc( ~zidx )         = num( ~zidx ) ./ den( ~zidx );
   
end

return   %calc_pearson_local
%------------------------------------------------------------------------%

