% LtL_make_Figure3           Generate figure 3
%
% Call: 
%               LtL_make_Figure3( tablename )
% Gets:
%               tablename   path the table levi2024_dataset
%
% 
% Calls:
%               LtL_calc_mr_svr
%
% Reference:
%               Levi, Aviv, and Stark, 2024 
%               "Learning to learn: Single session acquisition of new rules by freely moving mice", 
%               PNAS Nexus

% 27-Mar-24 ES and AL

% last update
% 26-may-24


function  LtL_make_Figure3( tablename )

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
mouse_num =  table2array(AA(:,3));
block_ordinare =  table2array(AA(:,4));
rule_num =  table2array(AA(:,6));
session_ordiante = table2array(AA(:,8));
SSL = table2array(AA(:,21));
MSL = table2array(AA(:,22));
intra_distance = table2array(AA(:,9));
inter_dist = table2array(AA(:,10));
for i = 1:length(inter_dist)
    temp = inter_dist{i};
    temp = str2double(temp);
    inter_distance(i) =  temp;
end
inter_distance = inter_distance';

% plot figure 3A - Success rates for individual rules:
urule = unique(rule_num);
nrule = length(urule);
rule_succ = NaN(300,6);
inter_rule = NaN(300,6);
intra_rule = NaN(300,6);
experience_rule = NaN(300,6);
SSL_rule = NaN(300,6);
MSL_rule = NaN(300,6);
mnum = NaN(300,6);
rule_nums = NaN(300,6);

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
    succ_per_block_use = [];
    inter_per_block_use = [];
    intra_per_block_use = [];
    SSL_per_block_use = [];
    MSL_per_block_use = [];
    exp_per_block_use = [];
    mnum_use = [];
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
end

figure;
boxplot(rule_succ, 'notch', 'on')
set(gca, 'TickDir','out','Box','off')
axis square
ylim([0 1.1])
ylabel('Success rate')
alines_local(0.5,'y','LineStyle', '--');
title('Figure 3A')

[ ~, ~, stats ] = kruskalwallis(rule_succ);
AA =  multcompare( stats );
pv = AA(:,6);

% plot figure 3B - Success rates vs. animal experience,
% quantified by the session ordinate
all_exp = experience_rule(:);
all_succ = rule_succ(:);
nanidx = isnan(all_exp);
all_exp = all_exp(~nanidx);
all_succ = all_succ(~nanidx);

uexp = unique(all_exp(~isnan(all_exp)));
nexp = length(uexp);
mean_succ = NaN(nexp,1);
sem_succ = NaN(nexp,1);

for i = 1:nexp
    idx = ismember(all_exp,uexp(i));
    mean_succ(i) = mean(all_succ(idx));
    sem_succ(i) = calc_sem_local(all_succ(idx));
end

figure;
scatter(uexp,mean_succ,'.')
lsline;
hold on
errorbar(uexp,mean_succ,sem_succ, 'o')
set(gca, 'TickDir','out','Box','off')
axis square
ylim([0.4 1])
ylabel('Success rate')
xlabel('Session ordinate')
alines_local(0.5,'y','LineStyle', '--');
title('Figure 3B1')

[cc, pv] = calc_spearman_local(all_exp,all_succ,1000);

%And only for SSL and MSL:
all_m = mnum(:);
all_rule = rule_nums(:);
all_SSL = SSL_rule(:);
all_MSL = MSL_rule(:);
all_SSL = all_SSL(~nanidx);
all_MSL = all_MSL(~nanidx);
all_rule = all_rule(~nanidx);
all_m=all_m(~nanidx);
%keep only first session of each MSL
mrule = [all_rule all_m];
umrule = unique(mrule, 'rows');
first_MSL = false(length(all_MSL),1);
for i = 1:length(umrule)
    idx = find(sum( mrule == umrule(i,:),2) == 2 & all_MSL==1);
    minsess = min(all_exp(idx));
    idx2 = idx(all_exp(idx)==minsess);
    first_MSL(idx2) = true;
end

idxSL = all_SSL|first_MSL;
expSL = all_exp(idxSL);
succSL = all_succ(idxSL);
uexp = unique(expSL);
nexp = length(uexp);
mean_succ = NaN(nexp,1);
sem_succ = NaN(nexp,1);
isSSL = all_SSL(idxSL);

for i = 1:nexp
    idx = ismember(expSL,uexp(i));
    mean_succ(i) = mean(succSL(idx));
    sem_succ(i) = calc_sem_local(succSL(idx));
end

figure;
scatter(uexp,mean_succ,'.')
lsline;
hold on
errorbar(uexp,mean_succ,sem_succ, 'o')
set(gca, 'TickDir','out','Box','off')
axis square
ylim([0.4 1])
ylabel('Success rate')
xlabel('Session ordinate')
alines_local(0.5,'y','LineStyle', '--');
title('Figure 3B2')
[cc, pv] = calc_spearman_local(expSL,succSL,1000);


% plot figure 3C - Success rates vs. DC difficulty, 
% quantified by the intra-criterion distance metric 
all_intra = intra_rule(:);
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
lsline;
hold on
errorbar(uintra,mean_succ,sem_succ, 'o')
set(gca, 'TickDir','out','Box','off')
axis square
ylim([0.4 1])
ylabel('Success rate')
xlabel('Intra-criterion distance')
alines_local(0.5,'y','LineStyle', '--');
title('Figure 3C1')
[cc, pv] = calc_spearman_local(all_intra,all_succ,1000);

%And only for SSL and MSL:

intraSL = all_intra(idxSL);
uintra = unique(intraSL);
nintra = length(uintra);
mean_succ = NaN(nintra,1);
sem_succ = NaN(nintra,1);

for i = 1:nintra
    idx = ismember(intraSL,uintra(i));
    mean_succ(i) = mean(succSL(idx));
    sem_succ(i) = calc_sem_local(succSL(idx));
end

figure;
scatter(uintra,mean_succ,'.')
hold on
errorbar(uintra,mean_succ,sem_succ, 'o')
set(gca, 'TickDir','out','Box','off')
axis square
ylim([0.4 1])
xlim([0 0.7])
lsline;

ylabel('Success rate')
xlabel('Session ordinate')
alines_local(0.5,'y','LineStyle', '--');
title('Figure 3C2')
[cc, pv] = calc_spearman_local(intraSL,succSL,1000);



% plot figure 3D - Success rates vs. the similarity
% between the new and the previous criterion
all_inter = inter_rule(:);
all_inter = all_inter(~nanidx);
all_inter = round(all_inter,5);
uinter = unique(all_inter(~isnan(all_inter)));
ninter = length(uinter);
mean_succ = NaN(ninter,1);
sem_succ = NaN(ninter,1);
for i = 1:ninter
    idx = ismember(all_inter,uinter(i));
    mean_succ(i) = mean(all_succ(idx));
    sem_succ(i) = calc_sem_local(all_succ(idx));
end

figure;
scatter(uinter,mean_succ,'.')
ylim([0.4 1])
xlim([0 0.9])
lsline;
hold on
errorbar(uinter,mean_succ,sem_succ, 'o')
set(gca, 'TickDir','out','Box','off')
axis square
ylim([0.4 1])
xlim([0 0.9])
ylabel('Success rate')
xlabel('Inter-criterion distance')
alines_local(0.5,'y','LineStyle', '--');
title('Figure 3D1')
[cc, pv] = calc_spearman_local(all_inter,all_succ,1000);

%And only for SSL and MSL:

interSL = all_inter(idxSL);
uinter = unique(interSL);
ninter = length(uinter);
mean_succ = NaN(ninter,1);
sem_succ = NaN(ninter,1);

for i = 1:ninter
    idx = ismember(interSL,uinter(i));
    mean_succ(i) = mean(succSL(idx));
    sem_succ(i) = calc_sem_local(succSL(idx));
end

figure;
scatter(uinter,mean_succ,'.')
ylim([0.4 1])
xlim([0 0.9])
lsline;
hold on
errorbar(uinter,mean_succ,sem_succ, 'o')
set(gca, 'TickDir','out','Box','off')
axis square
ylabel('Success rate')
xlabel('Session ordinate')
alines_local(0.5,'y','LineStyle', '--');
title('Figure 3D2')
[cc, pv] = calc_spearman_local(interSL,succSL,1000);


% plot figure 3E - cc-s and partial rank correlation coefficients (pcc) 
% between success rate and the features described in 3B-3D
Y = [expSL intraSL interSL];
figure;
[ pcc, cc, ~, ~, R2 ] = calc_mr_local( succSL,Y);
barwerror_local(1:3,cc(:,1));
hold on
h= barwerror_local(1:3,pcc(:,1),pcc(:,2),[0 0 0]);
set(h,'FaceColor','none')
set(h,'EdgeColor','k')
set(h,'LineWidth' ,2)
title(sprintf('Figure 3E \n R^2 = %0.2g\n pv = %0.2g, %0.2g, %0.2g',R2(1),pcc(:,3)));
axis square
set( gca, 'tickdir', 'out', 'box', 'off' )    
ylim([-0.5 0.5])
ylabel('Rank correlation coefficient')
set(gca,'XTickLabel',{'session ordinate'; 'intra-criterion distance'; 'inter-criterion distance'})


% plot figure 3F - Variance in block success rate (R2) explained by
% cross-validated support vector regression models
Y = [intraSL interSL expSL];
mat = [succSL Y];
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
    [ R2, res, fig ]            =  LtL_calc_mr_svr( mat( :, 2 : 4 ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );
    R2_svr( i )                 = R2;
end
        
% model comparison (importance):
R2_svr_models                   = NaN( nreps, 3 );
for i                           = 1 : nreps
    fprintf( 1, '%d ', i )
    if mod( i, 10 ) == 0
        fprintf( 1, '\n' )
    end
    [ R2_1, res, fig ]          =  LtL_calc_mr_svr( mat( :, [ 2 3 ] ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );
    [ R2_2, res, fig ]          =  LtL_calc_mr_svr( mat( :, [ 3 4 ] ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );
    [ R2_3, res, fig ]          =  LtL_calc_mr_svr( mat( :, [ 2 4 ] ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );
    R2_svr_models( i, : )       = [ R2_1 R2_2 R2_3];
end

figure
boxplot( R2_svr_models, 'notch', 'on' )
set( gca, 'tickdir', 'out', 'box', 'off' )
alines_local(median(R2_svr),'y');
axis square
ylim([0 0.5])
pv = [];
for i = 1:3
    pv(i) = utest_local( R2_svr_models( :, i ), R2_svr );
end
title(sprintf('Figure 3F \n SVR, pv = %0.2g,%0.2g,%0.2g',pv));
ylabel('Success rate R^2')
set(gca,'XTickLabel',{'session ordinate removed';'Intra-criterion removed';'inter-criterion removerd'});

% plot figure 3G -Accuracy in predicting SSL using 
% cross-validated support vector classification.
mat = [isSSL Y];
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
    [ R2, res, fig ]            =  LtL_calc_mr_svr( mat( :, 2 : 4 ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );
    R2_svr( i )                 = R2;
end
        
% model comparison (importance):
R2_svr_models                   = NaN( nreps, 3 );
for i                           = 1 : nreps
    fprintf( 1, '%d ', i )
    if mod( i, 10 ) == 0
        fprintf( 1, '\n' )
    end
    [ R2_1, res, fig ]          =  LtL_calc_mr_svr( mat( :, [ 2 3 ] ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );
    [ R2_2, res, fig ]          =  LtL_calc_mr_svr( mat( :, [ 3 4 ] ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );
    [ R2_3, res, fig ]          =  LtL_calc_mr_svr( mat( :, [ 2 4 ] ), mat( :, 1 ), 'nfold', 10, 'model', sv_model, 'preprocess', sv_preprocess );
    R2_svr_models( i, : )       = [ R2_1 R2_2 R2_3 ];
end


figure
boxplot( R2_svr_models, 'notch', 'on' )
set( gca, 'tickdir', 'out', 'box', 'off' )
alines_local(median(R2_svr),'y','LineStyle','--');
axis square
ylim([0.5 1.1])
pv = [];
for i = 1:3
    pv(i) = utest_local( R2_svr_models( :, i ), R2_svr );
end
title(sprintf('Figure 3G \n SVM, pv = %0.2g,%0.2g,%0.2g',pv));
ylabel('SSL accuracy')
set(gca,'XTickLabel',{'session ordinate removed';'Intra-criterion removed';'inter-criterion removerd'});

return           %LtL_make_Figure3

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
