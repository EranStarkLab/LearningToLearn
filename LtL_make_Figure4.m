% LtL_make_Figure4           Generate figure 4 
%
% Call: 
%               LtL_make_Figure4( tablename )
% Gets:
%               tablename   path the table levi2024_dataset
%
% 
% Calls:
%               nothing
%
% Reference:
%               Levi, Aviv, and Stark, 2024 
%               "Learning to learn: Single session acquisition of new rules by freely moving mice", 
%               PNAS Nexus

% 28-Apr-24 ES and AL

% last update
% 25-may-24

function  LtL_make_Figure4( tablename )

% argument handling
nargs                           = nargin;
if nargs < 1 || isempty( tablename )
    error( 'Path to dataset is required' )
end

% load sheet2 from the table:
AA = readtable(tablename,'Sheet',2,'VariableNamingRule','preserve');

if isempty(AA)
    error('Dataset is not available' )
end

correct_trial = table2array(AA(:,2));
mouse_num =  table2array(AA(:,3));
block_ordinare =  table2array(AA(:,4));
rule_num =  table2array(AA(:,6));
session_ordiante = table2array(AA(:,8));
SSL = table2array(AA(:,11));
MSL = table2array(AA(:,12));
intra_distance = table2array(AA(:,9));
inter_dist = table2array(AA(:,10));
inter_distance = [];
for i = 1:length(inter_dist)
    temp = inter_dist{i};
    temp = str2double(temp);
    inter_distance(i) =  temp;
end
inter_distance = inter_distance';

urule = unique(rule_num);
nrule = length(urule);
rule_succ = NaN(300,6);
inter_rule = NaN(300,6);
intra_rule = NaN(300,6);
experience_rule = NaN(300,6);
SSL_rule = NaN(300,6);
MSL_rule = NaN(300,6);
mnum = NaN(300,6);
rnum = NaN(300,6);

for i = 1:nrule
    idx = ismember(rule_num,urule(i));
    block_vec = block_ordinare(idx);
    sess_vec = session_ordiante(idx);
    usess = unique(sess_vec);
    succ_per_block = correct_trial(idx);
    inter_per_block = inter_distance(idx);
    intra_per_block = intra_distance(idx);
    SSL_per_block = SSL(idx);
    MSL_per_block = MSL(idx);
    exp_per_block = session_ordiante(idx);
    mnum_per_block = mouse_num(idx);
    rnum_per_block = rule_num(idx);
    succ_per_block_use = [];
    inter_per_block_use = [];
    intra_per_block_use = [];
    SSL_per_block_use = [];
    MSL_per_block_use = [];
    exp_per_block_use = [];
    mnum_use = [];
    rnum_use = [];
    l=0;
    for j = 1:length(usess)
        idx2 = ismember(sess_vec,usess(j));
        block_vec2 = block_vec(idx2);
        fi = find( diff( block_vec2 ) );
        ei = [ fi; length( block_vec2 ) ];
        si = [ 1; fi + 1 ];
        mat  = [ si ei ];
        nblocks = size( mat, 1 );
        succ_per_block2 = succ_per_block(idx2);
        inter_per_block2 = inter_per_block(idx2);
        intra_per_block2 = intra_per_block(idx2);
        SSL_per_block2 = SSL_per_block(idx2);
        MSL_per_block2 = MSL_per_block(idx2);
        exp_per_block2 = exp_per_block(idx2);
        mnum_per_block2 = mnum_per_block(idx2);
        rnum_per_block2 = rnum_per_block(idx2);
        for k = 1:nblocks
            l = l+1;            
            idx3 = mat( k, 1 ) : mat( k, 2 );
            succ_per_block_use(l) = mean(succ_per_block2(idx3));
            inter_per_block_use(l) = mean(inter_per_block2(idx3));
            intra_per_block_use(l) = mean(intra_per_block2(idx3));
            SSL_per_block_use(l) = mean(SSL_per_block2(idx3));
            MSL_per_block_use(l) = mean(MSL_per_block2(idx3));
            exp_per_block_use(l) = mean(exp_per_block2(idx3));
            mnum_use(l) = mean(mnum_per_block2(idx3));
            rnum_use(l) = mean(rnum_per_block2(idx3));
        end
    end
    numblock = length(succ_per_block_use);
    rule_succ(1:numblock,i) = succ_per_block_use;
    rule_nums(1:numblock,i) = i;
    inter_rule(1:numblock,i) = inter_per_block_use;
    intra_rule(1:numblock,i) = intra_per_block_use;
    experience_rule(1:numblock,i) = exp_per_block_use;
    SSL_rule(1:numblock,i) = SSL_per_block_use;
    MSL_rule(1:numblock,i) = MSL_per_block_use;
    mnum(1:numblock,i) = mnum_use;
    rnum(1:numblock,i) = rnum_use;
end

% plot figure 4C - Success rates for individual DCs.
% for rule#7, keep only the last session for each mouse:
rule7 =  NaN(300,1);
rule7(1:15,1) = rule_succ([56:66 100:103],1);
rule_succ2 = [rule7 rule_succ(:,2:5)];
figure;
boxplot(rule_succ2, 'notch', 'on')
set(gca, 'TickDir','out','Box','off')
axis square
ylim([0 1.1])
ylabel('Success rate')
alines_local(0.5,'y','LineStyle', '--');
title('Figure 4C')

[ ~, ~, stats ] = kruskalwallis(rule_succ);
AA =  multcompare( stats );
pv = AA(:,6);

% plot figure 4D - Success rate vs. DC difficulty.
all_succ = rule_succ2(:);
all_intra = intra_rule(:);
nanidx = isnan(all_succ);
all_succ = all_succ(~nanidx);
all_intra = all_intra(~nanidx);
all_intra = round(all_intra,5);
uintra = unique(all_intra);
nintra = length(uintra);
mean_succ = NaN(nintra,1);
sem_succ = NaN(nintra,1);
for i = 1:nintra
    idx = ismember(all_intra,uintra(i));
    mean_succ(i) = mean(all_succ(idx));
    sem_succ(i) = calc_sem_local(all_succ(idx));
end

figure;
scatter(uintra,mean_succ,'.')
ylim([0.4 1])
xlim([0 0.7])
lsline;
hold on
errorbar(uintra,mean_succ,sem_succ, 'o')
set(gca, 'TickDir','out','Box','off')
axis square
ylabel('Success rate')
xlabel('Intra-criterion distance')
alines_local(0.5,'y','LineStyle', '--');
title('Figure 4D')
text(uintra,mean_succ+0.1,{'11';'10';'9';'8';'7'})

[cc, pv] = calc_spearman_local(all_intra,all_succ,1000);

% plot figure 4E - Success rate vs. the
% similarity between the new and the previous DC.
all_succ = rule_succ2(:);
all_num = rnum(:);
all_inter = inter_rule(:);
nanidx = isnan(all_inter);
all_succ = all_succ(~nanidx);
all_inter = all_inter(~nanidx);
all_num = all_num(~nanidx);
all_inter = round(all_inter,5);
unum = unique(all_num(~isnan(all_num)));
nnum = length(unum);
mean_succ = NaN(nnum,1);
sem_succ = NaN(nnum,1);
uinter = all_inter(1:nnum,1);
for i = 1:nnum
    idx = ismember(all_num,unum(i));
    mean_succ(i) = mean(all_succ(idx));
    sem_succ(i) = calc_sem_local(all_succ(idx));
end

figure;
scatter(uinter,mean_succ,'.')
ylim([0.4 1])
xlim([0 0.9])
hold on
errorbar(uinter,mean_succ,sem_succ, 'o')
set(gca, 'TickDir','out','Box','off')
axis square
ylabel('Success rate')
xlabel('Inter-criterion distance')
alines_local(0.5,'y','LineStyle', '--');
title('Figure 4E')

return      %LtL_make_Figure4

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
