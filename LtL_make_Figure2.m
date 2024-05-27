% LtL_make_Figure2           Generate figure 2 
%
% Call: 
%               LtL_make_Figure2( tablename )
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

% 17-Mar-24  ES and AL

% last update
% 25-may-24

function LtL_make_Figure2( tablename )

% argument handling
nargs                           = nargin;
if nargs < 1 || isempty( tablename )
    error( 'Path to dataset is required' )
end

% load sheet1 from the table:
AA = readtable(tablename,'Sheet',1,'VariableNamingRule','preserve');

if isempty(AA)
    error('Dataset is not available' )
end

correct_trial = table2array(AA(:,2));
mouse_num =  table2array(AA(:,3));
block_ordinare =  table2array(AA(:,4));
rule_num =  table2array(AA(:,6));
rule_ordiante= table2array(AA(:,7));
session_ordiante = table2array(AA(:,8));
SSL = table2array(AA(:,21));
MSL = table2array(AA(:,22));

% plot figure 2A - example of SSL session from m4:
idxm = mouse_num == 4;
idxr = rule_num == 6;
idxs = session_ordiante == 11;

succ = correct_trial(idxm&idxr&idxs);
block_vec = block_ordinare(idxm&idxr&idxs);
numblcok = max(block_vec);
mean_succ = NaN(numblcok,1);
sem_succ = NaN(numblcok,1);

for i = 1:numblcok
    idx = block_vec == i;
    mean_succ(i) = mean(succ(idx));
    sem_succ(i) = calc_sem_local(succ(idx));
end

%fit linear line
P1 = polyfit(1:numblcok,mean_succ,1);
x = 0:0.01:(numblcok+1);
y = P1(2) + P1(1)*x;

%plot
figure;
errorbar(1:numblcok,mean_succ,sem_succ)
hold on
plot(x,y,'--k')
set(gca, 'TickDir','out','Box','off')
axis square
ylim([0 1])
xlim([0 numblcok+1])
xlabel('Block number')
ylabel('Success rate')
alines_local(0.5,'y','LineStyle', '--');
title('Figure 2A')

% plot figure 2B - example of MSL session from m4:
idxm = mouse_num == 4;
idxr = rule_num == 2;
idxs1 = session_ordiante == 6; 
idxs2 = session_ordiante==7;

succ = correct_trial(idxm&idxr&(idxs1|idxs2));
block_vec1 = block_ordinare(idxm&idxr&idxs1);
block_vec2 = block_ordinare(idxm&idxr&idxs2);
numblcok = max(block_vec1);
block_vec2 = block_vec2+numblcok;
block_vec = [block_vec1;block_vec2];
numblcok = max(block_vec);
mean_succ = NaN(numblcok,1);
sem_succ = NaN(numblcok,1);

for i = 1:numblcok
    idx = block_vec == i;
    mean_succ(i) = mean(succ(idx));
    sem_succ(i) = calc_sem_local(succ(idx));
end

%fit linear line
P1 = polyfit(1:numblcok,mean_succ,1);
x = 0:0.01:(numblcok+1);
y = P1(2) + P1(1)*x;

%plot
figure;
errorbar(1:numblcok,mean_succ,sem_succ)
hold on
plot(x,y,'--k')
set(gca, 'TickDir','out','Box','off')
axis square
ylim([0 1])
xlim([0 numblcok+1])
xlabel('Block number')
ylabel('Success rate')
alines_local(0.5,'y','LineStyle', '--');
title('Figure 2B')

% plot figure 2C - Calculate slopes of SSL and MSL
umouse = unique(mouse_num);
nummouse = length(umouse);
slope_SSL = [];
slope_MSL = [];
Intercepts_SSL = [];
Intercepts_MSL = [];
SSL_sessions = {};
l=0;
%For SSL:
for i = 1:nummouse
    idx = mouse_num == umouse(i);
    succ_SSL = correct_trial(idx&SSL);
    block_SSL = block_ordinare(idx&SSL);
    session_ord_SSL = session_ordiante(idx&SSL);
    usess = unique(session_ord_SSL);
    nsess = length(usess);
    for j = 1:nsess
        idx = session_ord_SSL==usess(j);
        if sum(idx) == 0
            continue
        end
        blockj = block_SSL(idx);
        succ_SSLj = succ_SSL(idx);
        numblock = max(blockj);
        succ_per_block = [];
        for k = 1:numblock
            succ_per_block(k) = mean(succ_SSLj(blockj==k));
        end
        l=l+1;
        P1 = polyfit(1:numblock,succ_per_block,1);
        slope_SSL(l) = P1(1);
        Intercepts_SSL(l) = P1(2);
        SSL_sessions{l} = succ_per_block;
    end
end

figure;
title('Figure 2C')
subplot(2,2,2)
histogram(slope_SSL,5,'DisplayStyle','stairs')
set(gca, 'TickDir','out','Box','off')
axis square
xlim([-0.07 0.07])
ylabel('Rules')
xlabel('Slope')
alines_local(0,'x','LineStyle','--','Color',[0 0 0]);
alines_local(mean(slope_SSL),'x','LineStyle','--','Color',[0 0 1]);
pv = signrank(slope_SSL);
title(sprintf('pv = %0.2g',pv));

subplot(2,2,4)
histogram(Intercepts_SSL,5,'DisplayStyle','stairs')
set(gca, 'TickDir','out','Box','off')
axis square
xlim([0 1])
ylabel('Rules')
xlabel('Initial success rate')
alines_local(0.5,'x','LineStyle','--','Color',[0 0 0]);
alines_local(mean(Intercepts_SSL),'x','LineStyle','--','Color',[0 0 1]);
pv = signrank(Intercepts_SSL,0.5);
title(sprintf('pv = %0.2g',pv));

%For MSL:
MSL_sessions = {};
l=0;
for i = 1:nummouse
    idx = mouse_num == umouse(i);
    succ_MSL = correct_trial(idx&MSL);
    block_MSL = block_ordinare(idx&MSL);
    session_ord_MSL = session_ordiante(idx&MSL);
    rule_ord_MSL = rule_ordiante(idx&MSL);
    urule = unique(rule_ord_MSL);
    nrule = length(urule);
    for j = 1:nrule
        idx = rule_ord_MSL==urule(j);
        if sum(idx) == 0
            continue
        end
        session_ord_MSLj = session_ord_MSL(idx);
        usess_in_MSL = unique(session_ord_MSLj);
        nsess_in_MSL = length(usess_in_MSL);
        blockj = block_MSL(idx); 
        blockj_all = [];
        max_block = 0;
        for n = 1:nsess_in_MSL
            idx2 = session_ord_MSLj == usess_in_MSL(n);
            blocks_temp = blockj(idx2);
            blockj_all = [blockj_all;blocks_temp+max_block];
            max_block = max(blocks_temp);
        end
        succ_MSLj = succ_MSL(idx);
        numblock = max(blockj_all);
        succ_per_block = [];
        for k = 1:numblock
            succ_per_block(k) = mean(succ_MSLj(blockj_all==k));
        end
        l=l+1;
        P1 = polyfit(1:numblock,succ_per_block,1);
        slope_MSL(l) = P1(1);
        Intercepts_MSL(l) = P1(2);
        MSL_sessions{l} = succ_per_block;
    end
end
subplot(2,2,1)
histogram(slope_MSL,5,'DisplayStyle','stairs')
set(gca, 'TickDir','out','Box','off')
axis square
xlim([-0.07 0.07])
ylabel('Rules')
xlabel('Slope')
alines_local(0,'x','LineStyle','--','Color',[0 0 0]);
alines_local(mean(slope_MSL),'x','LineStyle','--','Color',[0 0 1]);
pv = signrank(slope_MSL);
title(sprintf('pv = %0.2g',pv));

subplot(2,2,3)
histogram(Intercepts_MSL,5,'DisplayStyle','stairs')
set(gca, 'TickDir','out','Box','off')
axis square
xlim([0 1])
ylabel('Rules')
xlabel('Initial success rate')
alines_local(0.5,'x','LineStyle','--','Color',[0 0 0]);
alines_local(mean(Intercepts_MSL),'x','LineStyle','--','Color',[0 0 1]);
pv = signrank(Intercepts_MSL,0.5);
title(sprintf('pv = %0.2g',pv));

% plot figure 2D -Time warp for SSL and MSL
numSSL = length(SSL_sessions);
numMSL = length(MSL_sessions);
SSLwarp = NaN(numSSL,10);
for i = 1:numSSL
    temp = SSL_sessions{i};
    Vout = interp1(1:1/(length(temp)-1):2,temp,1:1/9:2,'spline');
    SSLwarp(i,:) = Vout; 
end
MSLwarp = NaN(numMSL,10);
for i = 1:numMSL
    temp = MSL_sessions{i};
    Vout = interp1(1:1/(length(temp)-1):2,temp,1:1/9:2,'spline');
    MSLwarp(i,:) = Vout; 
end

SSLwarp_plot = mean(SSLwarp);
SSLwarp_sem = calc_sem_local(SSLwarp);
MSLwarp_plot = mean(MSLwarp);
MSLwarp_sem = calc_sem_local(MSLwarp);

figure;
errorbar(0.1:0.1:1,SSLwarp_plot,SSLwarp_sem)
hold on
errorbar(0.1:0.1:1,MSLwarp_plot,MSLwarp_sem)
xlim([0 1])
ylim([0.4 0.9])
set(gca, 'TickDir','out','Box','off')
axis square
alines_local(0.5,'y','LineStyle','--','Color',[0 0 0]);
ylabel('Success rate')
xlabel('Fraction of learning passed')
pv = NaN(10,1);
title('Figure 2D')
for i = 1:10
    pv(i) = utest_local(SSLwarp(:,i),MSLwarp(:,i));
end

% plot figure 2E -Success rates of the first one,
% two, three… or ten testing trials of SSL and MSL sessions
succ_rt = table2array(AA(:,11:20));
succ_rt_SSL = succ_rt(logical(SSL),:);
succ_rt_MSL = succ_rt(logical(MSL),:);
diff_rule_num_SSL = [0;diff(rule_num(logical(SSL),:))];
idx_diff_SSL = [1;find(~ismember(diff_rule_num_SSL,0))];
diff_rule_num_MSL = [0;diff(rule_num(logical(MSL),:))];
idx_diff_MSL = [1;find(~ismember(diff_rule_num_MSL,0))]; 

mean_SSL_tr = NaN(10,1);
sem_SSL_tr = NaN(10,1);
mean_MSL_tr = NaN(10,1);
sem_MSL_tr = NaN(10,1);
for i = 1:10
    mean_SSL_tr(i) = mean(succ_rt_SSL(idx_diff_SSL,i));
    sem_SSL_tr(i) = calc_sem_local(succ_rt_SSL(idx_diff_SSL,i));
    mean_MSL_tr(i) = mean(succ_rt_MSL(idx_diff_MSL,i));
    sem_MSL_tr(i) = calc_sem_local(succ_rt_MSL(idx_diff_MSL,i));    
end

figure;
errorbar(1:10,mean_SSL_tr,sem_SSL_tr,'o')
hold on
errorbar(1:10,mean_MSL_tr,sem_MSL_tr,'o')
xlim([0 11])
ylim([0 1])
set(gca, 'TickDir','out','Box','off')
axis square
alines_local(0.5,'y','LineStyle','--','Color',[0 0 0]);
ylabel('Success rate: first-block trials')
xlabel('Number of first-block trials')
title('Figure 2E')
pv = NaN(10,1);
for i = 1:10
    pv(i) = utest_local(succ_rt_SSL(idx_diff_SSL,i),succ_rt_MSL(idx_diff_MSL,i));
end

% plot figure 2F -Success rates in all non-first blocks of the first 
% sessions of newly encountered DCs vs. success rates in the 
% first 3 trials of the session

%for SSL:
idx_diff_SSL2 = [idx_diff_SSL; length(succ_rt_SSL)];

SSL_first3 = NaN(length(idx_diff_SSL2)-1,1);
SSL_rest_of_session = NaN(length(idx_diff_SSL2)-1,1);
SSL_rest_of_session_err = NaN(length(idx_diff_SSL2)-1,1);
SSL_trials = correct_trial(logical(SSL));

for i =1:length(idx_diff_SSL2)-1
    idx = idx_diff_SSL2(i):idx_diff_SSL2(i+1);
    SSL_first3(i) = mean(SSL_trials(idx(1):idx(3)));
    SSL_rest_of_session(i) =  mean(SSL_trials(idx(11):idx(end)));
    SSL_rest_of_session_err(i) = calc_sem_local(SSL_trials(idx(11):idx(end)));
end

%for MSL:
idx_diff_MSL2 = [idx_diff_MSL; length(succ_rt_MSL)];

MSL_first3 = NaN(length(idx_diff_MSL2)-1,1);
MSL_rest_of_session = NaN(length(idx_diff_MSL2)-1,1);
MSL_rest_of_session_err = NaN(length(idx_diff_MSL2)-1,1);
MSL_trials = correct_trial(logical(MSL));
for i =1:length(idx_diff_MSL2)-1
    idx = idx_diff_MSL2(i):idx_diff_MSL2(i+1);
    MSL_first3(i) = mean(MSL_trials(idx(1):idx(3)));
    MSL_rest_of_session(i) =  mean(MSL_trials(idx(11):idx(end)));
    MSL_rest_of_session_err(i) = calc_sem_local(MSL_trials(idx(11):idx(end)));
end

figure;
errorbar(SSL_first3,SSL_rest_of_session,SSL_rest_of_session_err,'o')
hold on
errorbar(MSL_first3,MSL_rest_of_session,MSL_rest_of_session_err,'o')
xlim([0 1])
ylim([0 1])
set(gca, 'TickDir','out','Box','off')
axis square
alines_local(0.5,'y','LineStyle','--','Color',[0 0 0]);
ylabel('Success rate: all non-first blocks')
xlabel('Success rate - first block, 3 trials')
[cc, pv] = calc_spearman_local([SSL_first3;MSL_first3],[SSL_rest_of_session;MSL_rest_of_session],1000);
title('Figure 2F')

% plot figure 2G -cc-s between success rate during the first few first-block 
% testing trials and all non-first blocks, plotted against the number of 
% first-block testing trials. 

%for SSL:
SSL_first = NaN(length(idx_diff_SSL2)-1,10);
for j=1:10
    for i =1:length(idx_diff_SSL2)-1
        idx = idx_diff_SSL2(i):idx_diff_SSL2(i+1);
        SSL_first(i,j) = mean(SSL_trials(idx(1):idx(j)));
    end
end
%for MSL:
MSL_first = NaN(length(idx_diff_MSL2)-1,10);
for j = 1:10
    for i =1:length(idx_diff_MSL2)-1
        idx = idx_diff_MSL2(i):idx_diff_MSL2(i+1);
        MSL_first(i,j) = mean(MSL_trials(idx(1):idx(j)));
    end
end

cc_all = NaN(10,1);
pv_all = NaN(10,1);
for j = 1:10
    [cc_all(j),pv_all(j)] = calc_spearman_local([SSL_first(:,j);MSL_first(:,j)],[SSL_rest_of_session;MSL_rest_of_session],1000);
end

figure;
scatter(1:10,cc_all,'.')
xlim([0 11])
ylim([0 0.7])
set(gca, 'TickDir','out','Box','off')
axis square
ylabel('cc')
xlabel('Number of first-block trials')
title('Figure 2G')

return        %LtL_make_Figure2

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

