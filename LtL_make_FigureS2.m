% LtL_make_FigureS2           Generate figure S2    
%
% Call: 
%               LtL_make_FigureS2( tablename )
% Gets:
%               tablename   path the table levi2024_dataset
%
% 
% Calls:
%               nothing
%
% Reference:
%               Levi, Aviv, Stark, 2024 
%               "Learning to learn: Single session acquisition of new rules by freely-moving mice", 
%               PNAS nexus

% 18-may-24 ES and AL

% last update
% 25-may-24

function  LtL_make_FigureS2( tablename )

% argument handling
nargs                           = nargin;
if nargs < 1 || isempty( tablename )
    error( 'Path to dataset is required' )
end

%load sheet1 from the table:
AA = readtable(tablename,'Sheet',1);

if isempty(AA)
    error('Dataset is not available' )
end

correct_trial = table2array(AA(:,2));
rule_num =  table2array(AA(:,6));
SSL = logical(table2array(AA(:,21)));
MSL = logical(table2array(AA(:,22)));
mouse_num =  table2array(AA(:,3));
block_ordinare =  table2array(AA(:,4));
session_ordiante = table2array(AA(:,8));
intra_distance = table2array(AA(:,9));
inter_dist = table2array(AA(:,10));
for i = 1:length(inter_dist)
    temp = inter_dist{i};
    temp = str2double(temp);
    inter_distance(i) =  temp;
end
inter_distance = inter_distance';
succ_rt = table2array(AA(:,11:20));

% plot figure S2A -Success rates
% during the first 1, 2, …, or 10 trials of SSL and MSL 
% sessions is correlated with intra-criterion distance.
% cc, rank correlation coefficient.

% find all the first sessions of MSL and SSL
succ_rt_use = [];
intra_to_use = [];
succ_per_block =[];
sess_ordin_per_block=[];
inter_per_block = [];
first_per_block = [];
isSSL_per_block = [];
umouse = unique(mouse_num);
nmouse = length(umouse);
l=0;
for i = 1:nmouse
    idx = ismember(mouse_num,umouse(i));
    SSL_idx = SSL(idx);
    MSL_idx = MSL(idx);
    session_idx = session_ordiante(idx);
    sess_SSL = unique(session_idx(SSL_idx));
    sess_MSL = unique(session_idx(MSL_idx));
    % keep only fist session for MSL:
    diffMSL = [2; diff(sess_MSL)];
    idxfirst = diffMSL>1;
    if ~isempty(sess_MSL)
    sess_MSL = sess_MSL(idxfirst);
    end
    if i == 6
        sess_MSL = [sess_MSL;13];
    end
    %find the success rate in the first trials and the intra-CD distance:
    all_first_sess = [sess_SSL ;sess_MSL];
    for j = 1:length(all_first_sess)
        idx2 = ismember(session_ordiante,all_first_sess(j));
        succ_rt_temp = succ_rt(idx&idx2,:);
        succ_rt_temp = unique (succ_rt_temp,'rows');
        succ_rt_use = [succ_rt_use;succ_rt_temp];
        intra_temp = intra_distance(idx&idx2);
        intra_temp = unique(intra_temp);
        intra_to_use = [intra_to_use;intra_temp];

        %find the session ordiante, intra-criterion
        % inter critertion and success rate per
        % block in the first session (use for fig S2B-F)
        blocks = block_ordinare(idx&idx2);
        ublocks = unique(blocks);
        nblocks = length(ublocks);
        for k = 1:nblocks
            l=l+1;
            idx3 = ismember(block_ordinare,ublocks(k));
            correct_trial_block = correct_trial(idx&idx2&idx3);
            succ_per_block(l) = mean(correct_trial_block);
            sess_ordin_per_block(l) = all_first_sess(j);
            intra_per_block(l) = intra_temp;
            inter_block = inter_distance(idx&idx2&idx3);
            inter_per_block(l) = mean(inter_block);
            first_three_block = succ_rt(idx&idx2&idx3,3);
            first_per_block(l) = mean(first_three_block);
            isSSL = SSL(idx&idx2&idx3);
            isSSL_per_block(l) = mean(isSSL);
        end
    end
end

figure; 
%calcularte the cc:
cc_all = NaN(10,1);
pv_all = NaN(10,1);
for j = 1:10
    [cc_all(j),pv_all(j)] = calc_spearman_local(succ_rt_use(:,j),intra_to_use,1000);
end

scatter(1:10,cc_all);
set(gca, 'TickDir','out','Box','off')
xlim([0 11]);
ylim([0 1]);
lsline
ylabel('cc between first-block success rate and intra-criterion distance')
xlabel('Number of first-block trials')
[cc,pv] = calc_spearman_local(1:10,cc_all,1000);
title(sprintf('figure S2A, cc=%0.2g, pv=%0.2g',cc,pv));

% plot figure S2B -Correlation matrix for session ordinate,
% intra- criterion distance, inter- criterion distance and the success rate during the
% first three trials.
Y = [sess_ordin_per_block' intra_per_block' inter_per_block' first_per_block'];
[C PP] = corrcoef(Y);
C = tril(C,-1);
C(C==0) = NaN;
C(logical(eye(size(C)))) = 1;
figure;
imagescnan(C)
set(gca,'YDir','Reverse')  
colormap('myjet')
set(gca, 'TickDir','out','Box','off');
set(gca,'YTick',[1,2,3,4]);
set(gca,'XTick',[1,2,3,4]);
set(gca,'XTickLabel',{'session ordinare';'intra-criterion';'inter-criterion';'first 3 trials'})
set(gca,'YTickLabel',{'session ordinare';'intra-criterion';'inter-criterion';'first 3 trials'})
colorbar
title('figure S2B')

% plot figure S2C -cc’s and pcc’s between success rate and session ordinate,
% intra- criterion distance, inter- criterion distance and the success rate during the
% first three trials.
Y = [sess_ordin_per_block' intra_per_block' inter_per_block' first_per_block'];
figure;
[ pcc, cc, ~, ~, R2 ] = calc_mr_local( succ_per_block',Y);
barwerror_local(1:4,cc(:,1));
hold on
h= barwerror_local(1:4,pcc(:,1),pcc(:,2),[0 0 0]);
set(h,'FaceColor','none')
set(h,'EdgeColor','k')
set(h,'LineWidth' ,2)
title(sprintf('Figure S2C \n R^2 = %0.2g\n pv = %0.2g, %0.2g, %0.2g, %0.2g',R2(1),pcc(:,3)));
axis square
set( gca, 'tickdir', 'out', 'box', 'off' )    
ylim([-0.5 0.5])
ylabel('Rank correlation coefficient')
set(gca,'XTickLabel',{'session ordinate'; 'intra-criterion distance'; 'inter-criterion distance';'first 3 trials'})


% plot figure S2D - Variance in block success rate (R^2) explained by cross-
%validated support vector regression models.
Y = [sess_ordin_per_block' intra_per_block' inter_per_block' first_per_block'];
mat = [succ_per_block' Y];
sv_preprocess                   = 'scale';
sv_model                        = 'svr';

% cross-validated SVR:
% get the variance for the R2's
nreps                           = 20;
R2_svr                          = NaN( nreps, 1 );
for i                           = 1 : nreps
    fprintf( 1, '%d ', i )
    if mod( i, 10 ) == 0
        fprintf( 1, '\n' )
    end
    [ R2, res, fig ]            =  LtL_calc_mr_svr( mat( :, 2 : 5 ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );
    R2_svr( i )                 = R2;
end
        
% model comparison (importance):
R2_svr_models                   = NaN( nreps, 4 );
for i                           = 1 : nreps
    fprintf( 1, '%d ', i )
    if mod( i, 10 ) == 0
        fprintf( 1, '\n' )
    end
    [ R2_1, res, fig ]          =  LtL_calc_mr_svr( mat( :, [ 3 4 5 ] ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );
    [ R2_2, res, fig ]          =  LtL_calc_mr_svr( mat( :, [ 2 4 5 ] ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );
    [ R2_3, res, fig ]          =  LtL_calc_mr_svr( mat( :, [ 2 3 5 ] ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );
    [ R2_4, res, fig ]          =  LtL_calc_mr_svr( mat( :, [ 2 3 4 ] ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );    
    R2_svr_models( i, : )       = [ R2_1 R2_2 R2_3 R2_4];
end

figure
boxplot( R2_svr_models, 'notch', 'on' )
set( gca, 'tickdir', 'out', 'box', 'off' )
alines_local(median(R2_svr),'y');
axis square
ylim([0 0.5])
pv = [];
for i = 1:4
    pv(i) = utest_local( R2_svr_models( :, i ), R2_svr );
end
title(sprintf('Figure S2D \n SVR, pv = %0.2g,%0.2g,%0.2g,%0.2g',pv));
ylabel('Success rate R^2')
set(gca,'XTickLabel',{'session ordinate removed';'Intra-criterion removed';'inter-criterion removerd';'first 3 trials removed'});

% plot figure S2E -Accuracy in predicting SSL using 
% cross-validated support vector classification.
mat = [isSSL_per_block' Y];
sv_preprocess                   = 'scale';
sv_model                        = 'svm';


% cross-validated SVM:
% get the variance for the R2's
nreps                           = 20;
R2_svr                          = NaN( nreps, 1 );
for i                           = 1 : nreps
    fprintf( 1, '%d ', i )
    if mod( i, 10 ) == 0
        fprintf( 1, '\n' )
    end
    [ R2, res, fig ]            =  LtL_calc_mr_svr( mat( :, 2 : 5 ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );
    R2_svr( i )                 = R2;
end
        
% model comparison (importance):
R2_svr_models                   = NaN( nreps, 4 );
for i                           = 1 : nreps
    fprintf( 1, '%d ', i )
    if mod( i, 10 ) == 0
        fprintf( 1, '\n' )
    end
    [ R2_1, res, fig ]          =  LtL_calc_mr_svr( mat( :, [ 3 4 5 ] ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );
    [ R2_2, res, fig ]          =  LtL_calc_mr_svr( mat( :, [ 2 4 5 ] ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );
    [ R2_3, res, fig ]          =  LtL_calc_mr_svr( mat( :, [ 2 3 5 ] ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );
    [ R2_4, res, fig ]          =  LtL_calc_mr_svr( mat( :, [ 2 3 4 ] ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );   
    R2_svr_models( i, : )       = [ R2_1 R2_2 R2_3 R2_4];
end


figure
boxplot( R2_svr_models, 'notch', 'on' )
set( gca, 'tickdir', 'out', 'box', 'off' )
alines_local(median(R2_svr),'y','LineStyle','--');
axis square
ylim([0.5 1.1])
pv = [];
for i = 1:4
    pv(i) = utest_local( R2_svr_models( :, i ), R2_svr );
end
title(sprintf('Figure S2E \n SVM, pv = %0.2g,%0.2g,%0.2g,%0.2g',pv));
ylabel('SSL accuracy')
set(gca,'XTickLabel',{'session ordinate removed';'Intra-criterion removed';'inter-criterion removerd';'first 3 trials removed'});

% plot figure S2F -Accuracy in predicting SSL using single-
%feature models.
% model comparison (importance):
R2_svr_models                   = NaN( nreps, 4 );
for i                           = 1 : nreps
    fprintf( 1, '%d ', i )
    if mod( i, 10 ) == 0
        fprintf( 1, '\n' )
    end
    [ R2_1, res, fig ]          =  LtL_calc_mr_svr( mat( :, [ 2 ] ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );
    [ R2_2, res, fig ]          =  LtL_calc_mr_svr( mat( :, [ 3 ] ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );
    [ R2_3, res, fig ]          =  LtL_calc_mr_svr( mat( :, [ 4 ] ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );
    [ R2_4, res, fig ]          =  LtL_calc_mr_svr( mat( :, [ 5 ] ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );   
    R2_svr_models( i, : )       = [ R2_1 R2_2 R2_3 R2_4];
end

figure
boxplot( R2_svr_models, 'notch', 'on' )
set( gca, 'tickdir', 'out', 'box', 'off' )
alines_local(median(R2_svr),'y','LineStyle','--');
axis square
ylim([0.5 1.1])
pv = [];
for i = 1:4
    pv(i) = utest_local( R2_svr_models( :, i ), R2_svr );
end
title(sprintf('Figure S2F \n SVM, pv = %0.2g,%0.2g,%0.2g,%0.2g',pv));
ylabel('SSL accuracy')
set(gca,'XTickLabel',{'only session ordinate';'only intra-criterion';'only inter-criterion';'only first 3 trials'});

return


%------------------------------------------------------------------------
%  y = calc_sem_local( x, dim )
% calculate the standard error of the mean of a vctor x
%------------------------------------------------------------------------%
% gets              X       data
%                   DIM     dimension {1} - sem of columns

function y = calc_sem_local( x, dim )

if isempty( x )
    y                           = NaN;
    return
end
if all( isnan( x( : ) ) )
    if nargin < 2 || isempty( dim )
        dim                     = 1;
    end
    y                           = nanmean( x, dim );
    return 
end 
if ~exist( 'dim', 'var' ) || isempty( dim )
    if any( size( x ) == 1 )
        x                       = x( : );
    end
    dim                         = 1;
end
y                               = nanstd( x, [], dim ) ./ sqrt( sum( ~isnan( x ), dim ) );
return    %calc_sem_local
%------------------------------------------------------------------------



%------------------------------------------------------------------------
% lh = alines_local( x, ax, varargin )
% plots seperators on current axis.
%------------------------------------------------------------------------%
% x         points to place separators
% ax        which axis to separate {'x'}
% the rest of the arguments are passed to line
%
% lh        handle to the lines
function lh = alines_local( x, ax, varargin )

if nargin < 2 || isempty( ax ), ax = 'x'; end
ax = lower( ax( 1 ) );

x = x( : );
nx = length( x );
ah = gca;
if ax == 'x'
    X = [ x x NaN * ones( nx, 1 ) ]';
    X = X( : );
    Y = repmat( [ get( ah, 'ylim' ) NaN ]', nx, 1 );
elseif ax == 'y'
    Y = [ x x NaN * ones( nx, 1) ]';
    Y = Y( : );
    X = repmat( [ get( ah, 'xlim' ) NaN ]', nx, 1 );
else
    return
end
lh = line( X, Y, varargin{:} );


return      %alines_local
%------------------------------------------------------------------------


%------------------------------------------------------------------------
% [ cc, pval ] = calc_spearman_local( x, y, nreps )
% calculate the spearman correlation
%------------------------------------------------------------------------%
% x       data matrix (2 columns)
% y       optional; then computes the CC between columns of x and y
% nreps   optional; then conducts a permutation test, nreps times
function [ cc, pval ] = calc_spearman_local( x, y, nreps )

% intialize output
cc                      = [];
pval                    = [];

% handle arguments
if nargin < 1 || isempty( x )
    return
end
if nargin < 2
    y                   = [];
end
if nargin < 3 || isempty( nreps )
    nreps               = 0;
end

% handle NaNs
nansX                   = isnan( x );
if isempty( y )
    nans                = sum( nansX, 2 ) > 0;
    x                   = x( ~nans, : );
elseif isequal( size( x ), size( y ) ) && size( x, 2 ) == 1
    nans                = sum( isnan( y ) | nansX, 2 ) > 0;
    x                   = x( ~nans, : );
    y                   = y( ~nans, : );
else
    % ignore NaNs for now, handled in calc_pearson_local
end

% rank-order and compute
x                       = rankcols_local( x );
y                       = rankcols_local( y );
cc                      = calc_pearson_local( x, y );

% prepare for significance testing
if nreps == 0 || isempty( x )
    if isempty( x ) && all( nans( : ) )
        cc              = NaN;
    end
    pval                = NaN;
    return
end
if isempty( y )                                 % simple case #1
    x1                  = x( :, 1 );
    x2                  = x( :, 2 );
elseif size( x, 2 ) == 1 && size( y, 2 ) == 1   % simple case #2
    x1                  = x;
    x2                  = y;
elseif size( x, 2 ) == size( y, 2 )             % recursion
    
    npairs              = size( x, 2 );
    cc                  = NaN * ones( 1, npairs );
    pval                = cc;
    for i               = 1 : npairs
        if all( isnan( x( :, i ) ) ) || all( isnan( y( :, i ) ) )
            cc( i )     = NaN;
            pval( i )   = NaN;
        else
            [ cc( i ), pval( i ) ] = calc_spearman_local( x( :, i ), y( :, i ), nreps );
        end
    end
    return
else                                            % aberrant matrices
    pval                = NaN * ones( size( cc ) );
    return
end

% estimate significance
v                       = ones( 1, nreps );
ccmix                   = calc_spearman_local( x1 * v, mixmat_local( x2 * v, 1 ) );
pval                    = ( sum( abs( cc ) <= abs( ccmix ) ) + 1 ) / ( nreps + 1 );

return %calc_spearman_local
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
% [ y, indx ]    = mixmat_local(x,dim,mode)
% mix matrix elements
%------------------------------------------------------------------------%
%   x       data matrix
%   matrix elements are mixed in dimension dim (1 = cols).
%   if is zero or left empty, all elements are mixed.
%   mode is 1 (w/o replacement) or 2 (with).
%   Y is same size as X.

function [ y, indx ]    = mixmat_local(x,dim,mode)

global resetClock
if isempty( resetClock )
    resetClock          = 1;
end
if resetClock
    rng( round( rand( 1 ) * sum( 100 * clock ) ), 'v4' )
    resetClock          = 0;
end

nargs                   = nargin;
if nargs < 1
    error( '1 argument')
end
if nargs < 2 || isempty(dim)
    dim                 = 1; 
end
if nargs < 3 || isempty(mode)
    mode                = 1; 
end

switch dim
    case 0
        isx             = size( x ); 
        x               = x( : ); 
        cdim            = 1;
    case 1
        cdim            = dim;
    case 2
        x               = x'; 
        cdim            = ( ~( dim - 1 ) ) + 1;
end

sx                      = size( x );
if ~prod( sx )
    y                   = zeros( sx );
    return
end
ridx                    = 0 : sx( 1 ) : ( sx( 2 ) - 1 ) * sx( 1 );

if mode == 1                                            % without replacement
    [ ~, p ]            = sort( rand( sx ), cdim );
elseif mode == 2                                        % with replacement
    p                   = ceil( rand( sx ) * sx( 1 ) );
end

indx                    = ridx(ones(sx(1),1),:) + p;

switch dim
    case 0
        y               = reshape( x( indx ), isx );
    case 1
        y               = x( indx );
    case 2
        y               = x( indx )';
end

return  %mixmat_local
%------------------------------------------------------------------------%


%------------------------------------------------------------------------
% [ pcc, cc, c, beta, R2 ] = calc_mr_local( y, x, rankf, nreps )
% calculate multiple regression
%------------------------------------------------------------------------%
% y           system output, n x 1
% x           system input, n x m
% rankf       {1}; flag for rank correlation
% nreps       {100}; number of bootstrap for SEM of contribution and R2
%
% returns               pcc         partial cc      m x 3: [ pcc    SD  p-values ]
%                       cc          cc              m x 3: [ cc     SD  p-values ]
%                       c           contribution    m x 2: [ c      SD ]
%                       beta        beta            m x 3: [ beta   SD  p-values ]
%                       R2          R2              1 x 3: [ R2     SD  p-value  ]
function [ pcc, cc, c, beta, R2 ] = calc_mr_local( y, x, rankf, nreps )

% arguments
nargs = nargin;
if nargs < 2 || isempty( y ) || isempty( x )
    error( 'missing inputs' )
end
[ n, m ]                        = size( x );
if n ~= size( y, 1 )
    error( 'input size mismatch' )
end
if nargs < 3 || isempty( rankf )
    rankf                       = 1;
end
if nargs < 4 || isempty( nreps )
    nreps                       = 100;
end

% rank
if rankf
    y                           = rankcols_local( y );
    x                           = rankcols_local( x );
end

% partial correlation, standardized regression 
[ pcc, beta ]                   = calc_pcc_local( y, x );

% correlation 
cc                            	= calc_pearson_local( y * ones( 1, m ), x )';

% significance (correlations, partial correlations)
pvals_cc                        = NaN( m, 1 );
pvals_pcc                       = NaN( m, 1 );
for i                           = 1 : m
    [ ~, pvals_cc( i) ]         = partialcorr( y, x( :, i ), ones( n, 1 ) );
    cidx                        = setdiff( 1 : m, i );
    [ ~, pvals_pcc( i ) ]       = partialcorr( y, x( :, i ), x( :, cidx ) );
end

% contribution and R2
c                               = beta .* cc; % C may be negative for some of the regressors
R2                              = sum( c );

% compute SEM for the contribution and R2
mix                             = mixmat_local( ( 1 : n )' * ones( 1, nreps ), 1, 2 );
cc_bs                           = NaN( m, nreps );
pcc_bs                          = NaN( m, nreps );
beta_bs                         = NaN( m, nreps );
for i                           = 1 : nreps 
    ridx                        = mix( :, i );
    [ pcc_bs( :, i ) ...
        , beta_bs( :, i ) ]     = calc_pcc_local( y( ridx, : ), x( ridx, : ) );
    cc_bs( :, i )               = calc_pearson_local( y( ridx, : ) * ones( 1, m ), x( ridx, : ) )';
end
sd_cc                           = calc_sem_local( cc_bs, 2 );
sd_pcc                          = calc_sem_local( pcc_bs, 2 );
sd_beta                         = calc_sem_local( beta_bs, 2 );
c_bs                            = pcc_bs .* beta_bs;
sd_c                            = calc_sem_local( c_bs, 2 );
sd_R2                           = calc_sem_local( sum( c_bs, 1 ) );

% p-value for full model
[ ~, ~, ~, ~, stts ]            = regress( y, [ ones( n, 1 ) x ] );

% organize output
cc                              = [ cc sd_cc pvals_cc ];
pcc                             = [ pcc sd_pcc pvals_pcc ];
beta                            = [ beta sd_beta pvals_pcc ];
c                               = [ c sd_c ];
R2                              = [ R2 sd_R2 stts( 3 ) ];

return          %calc_mr_local
%------------------------------------------------------------------------


%------------------------------------------------------------------------
% [ rp, beta ] = calc_pcc_local( y, X )
% calculate all maxmial order partial correlation coefficients and 
% standardized regression coefficients
%------------------------------------------------------------------------%
% y           system output, n x 1
% x           system input, n x m
function [ rp, beta ] = calc_pcc_local( y, X )

[ n, p ]       	= size( X );
y               = y - sum( y ) / n;
s               = sqrt( sum( y .* y, 1 ) / ( n - 1 ) );
if s == 0
    y           = zeros( n, 1 ); 
else
    y           = y / s; 
end
m               = sum( X, 1 ) / n;
d               = ones( n, 1 );
X               = X - d * m;
s               = sqrt( sum( X .* X, 1 ) / ( n - 1 ) );
zidx            = s == 0;
den             = d * s;
Z               = zeros( n, p );
Z( :, ~zidx )   = X( :, ~zidx ) ./ den( :, ~zidx );
I               = eye( p );
[ Q, R ]        = qr( Z, 0 );
beta            = ( R \ ( Q' * y ) );
Tstat           = beta ./ sqrt( diag( R \ I * ( R \ I )' ) ) / norm( y - Z * beta ) * sqrt( n - p );
rp              = Tstat ./ ( sqrt( Tstat .^ 2 + n - p ) );

return                  %calc_pcc_local
%------------------------------------------------------------------------


%------------------------------------------------------------------------
% [ bh, eh ] = barwerror_local( x, y, e, barcolor, barwidth, ebarcolor, ebarlinewidth, ebarwidth, orientation )
% plot a bar graph w/ error bars
%------------------------------------------------------------------------%
% x       (optional!)
% y       values
% e       error, the error bars will be 2e
function [ bh, eh ] = barwerror_local( x, y, e, barcolor, barwidth, ebarcolor, ebarlinewidth, ebarwidth, orientation )

% initialize output
bh                      = [];
eh                      = [];

% arguments
nargs                   = nargin;
if nargs < 2 || isempty( y )
    return
end
x                       = x( : );
y                       = y( : );
n                       = length( y );
if isempty( x )
    x                   = 1 : n;
end
if length( x ) ~= n
    return
end
if nargs < 3 || isempty( e )
    e                   = [];
end
e                       = e( : );
if ~isempty( e )
    if length( e ) ~= n
        return
    end
end
if nargs < 4 || isempty( barcolor ) || numel( barcolor ) ~= 3
    barcolor            = [ 0 0 1 ];
end
if nargs < 5 || isempty( barwidth )
    barwidth            = 0.8;
end
if nargs < 6 || isempty( ebarcolor )
    ebarcolor           = [ 0 0 0 ];
end
if nargs < 7 || isempty( ebarlinewidth )
    ebarlinewidth       = 1;
end
if nargs < 8 || isempty( ebarwidth )
    if n == 1
        ebarwidth       = barwidth / 2;
    else
        ebarwidth       = mean( diff( x ) ) * barwidth / 4;
    end
end
if nargs < 9 || isempty( orientation )
    orientation = 'vertical';
end
orientation             = lower( orientation( 1 ) );

% plot
if orientation == 'h'
    bh                  = barh( x, y, barwidth );
else
    bh                  = bar( x, y, barwidth );
end
set( bh, 'EdgeColor', barcolor, 'FaceColor', barcolor )
if length( x ) > 1
    if orientation == 'h'
        ylim( x( [ 1 end ] )' + ( barwidth + diff( x( 1 : 2 ) ) ) / 2 * [ -1 1 ] )
    else
        xlim( x( [ 1 end ] )' + ( barwidth + diff( x( 1 : 2 ) ) ) / 2 * [ -1 1 ] )
    end
end

if ~isempty( e )
    hstate = ishold;
    if ~hstate
        hold on
    end
    eh = make_cl_on_bar_local( x, y + e, y - e, ebarcolor, ebarwidth, ebarlinewidth, orientation );
    if ~hstate
        hold off
    end
end

return           %barwerror_local
%------------------------------------------------------------------------


%------------------------------------------------------------------------
% lh = make_cl_on_bar_local( x, ytop, ybot, color, width, linewidth, orientation )
%%------------------------------------------------------------------------%

function lh = make_cl_on_bar_local( x, ytop, ybot, color, width, linewidth, orientation )

% argument handling
nargs                   = nargin;
if nargs < 4 || isempty( color )
    color               = [ 0 0 1 ];
end
if nargs < 5 || isempty( width )
    width               = 0.4;
end
if nargs < 6 || isempty( linewidth )
    linewidth           = 1;
end
if nargs < 7 || isempty( orientation )
    orientation         = 'vertical';
end
orientation             = lower( orientation( 1 ) );

% prepare for plotting
x                       = x( : );
ytop                    = ytop( : );
ybot                    = ybot( : );
len                     = [ length( x ), length( ytop ), length( ybot ) ]; 
n                       = max( len );
if length( unique( length( len > 1 ) ) ) > 1
    error( 'mismatch' )
end
if len( 1 ) == 1
    x                   = repmat( x, [ n 1 ] );
end
if len( 2 ) == 1
    ytop                = repmat( ytop, [ n 1 ] );
end
if len( 3 ) == 1
    ybot                = repmat( ybot, [ n 1 ] );
end

% actually plot
lh                      = zeros( n, 1 );
for i                   = 1 : n
    if orientation == 'h'
        lh( i )         = line( [ ytop( i ) * [ 1 1 1 ] ybot( i ) * [ 1 1 1 ] ]...
            , [ x( i ) + [ -1 1 ] * width / 2 [ x( i ) x( i ) ] x( i ) + [ -1 1 ] * width / 2 ] );    
    else
        lh( i )         = line( [ x( i ) + [ -1 1 ] * width / 2 [ x( i ) x( i ) ] x( i ) + [ -1 1 ] * width / 2 ]...
            , [ ytop( i ) * [ 1 1 1 ] ybot( i ) * [ 1 1 1 ] ] );
    end
    set( lh( i ), 'color', color, 'linewidth', linewidth );
end

return              %make_cl_on_bar_local
%------------------------------------------------------------------------%

%------------------------------------------------------------------------
% [ p, U, method ] = utest_local( x1, x2, test )
% Mann-Whitney U-test
%------------------------------------------------------------------------%
% X1, X2      samples; may be unequal size
% TEST        test type:
%                               {0}     two sided, H1: X1 != X2
%                                1      one sided, H1: X1 >  X2
%                               -1      one sided, H1: X1 <  X2
function [ p, U, method ] = utest_local( x1, x2, test )
% handle arguments
nargs = nargin;
if nargs < 2
    if size( x1, 2 ) == 2
        x2 = x1( :, 2 );
        x1 = x1( :, 1 );
    else
        error( '2 arguments' )
    end
end
if nargs < 3 || isempty( test ), test = 0; end
% handle special cases
p = NaN;
U = NaN;
method = 'none'; 
if isequal( x1, x2 ), p = 1; return, end
if isempty( x1 ) || isempty( x2 ), p = 1; return, end

% sample sizes
x1 = x1( ~isnan( x1 ) );
x2 = x2( ~isnan( x2 ) );
x1 = x1( : );
x2 = x2( : );
n1 = length( x1 ); 
n2 = length( x2 );
n = n1 + n2;
if n1 == 0 || n2 == 0
    return
end
    
% Mann-Whitney U statistic

[ ranks, T ] = rankcols_local( [ x1; x2 ] );
R = sum( ranks( 1 : n1 ) );
C = n1 * n2 + n1 * ( n1 + 1 ) / 2 - R;              % p.429
U = max( C, n1 * n2 - C );

% significance

expcR = n1 * ( n + 1 ) / 2;
if ( n1 + n2 ) < 20
    if T                                            % exact binomial
        f = sum( nchoosek( ranks, n1 ), 2 );
        if R < expcR
            p = sum( f <= R ) / length( f );
        else 
            p = sum( f >= R ) / length( f );
        end
        method = 'exact';
    else                                            % use a table based on no-ties (perfect for that case only)
        eval( sprintf( 'load( ''utable.mat'', ''ptab_%d'' )', n ) )
        if R < expcR
            eval( sprintf( 'p = ptab_%d( ceil( R ), n1 );', n ) )
        else
            eval( sprintf( 'p = 1 - ptab_%d( ceil( R ) - 1, n1 );', n ) )
        end
        method = 'lookup';
    end
else                                                % approximate gaussian
    nom = U - n1 * n2 / 2;
    if T                                            % correction for ties, p.430
        den = n1 * n2 / 12 / n / ( n - 1 ) * ( n .^ 3 - n - T );
    else
        den = n1 * n2 / 12 * ( n + 1 );
    end
    ts = nom / sqrt( den );
    p = normcdf( -abs( ts ), 0, 1 );
    method = 'gaussian';
end

% translate p-value according to test type

if ~test
    p = 2 * p;
elseif ( test == 1 && R < expcR ) || ( test == -1 && R > expcR )
    p = 1 - p;
end

return  %utest_local
%------------------------------------------------------------------------%
