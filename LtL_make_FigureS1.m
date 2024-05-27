% LtL_make_FigureS1           Generate figure S1    
%
% Call: 
%               LtL_make_FigureS1( tablename )
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

% 18-may-24 ES and AL

% last update
% 25-may-24

function  LtL_make_FigureS1( tablename )

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
rule_ordiante= table2array(AA(:,7));
session_ordiante = table2array(AA(:,8));
SSL = logical(table2array(AA(:,21)));
MSL = logical(table2array(AA(:,22)));
first_session_learned = table2array(AA(:,23));
prob = 0.5;

% plot figure S1B -Mice learn a new DC
%  within a median of 90 testing trials. Data are from six 
% mice tested on 28 DCs.

umouse = unique(mouse_num);
nmouse = length(umouse);
numtr2succ=[];
mnum2succ =[];
issucc =[];
numtr = [];
nmou = [];
succperblock = [];
issuccperblock = [];
succperblockSSL = [];
isSSL = [];
succperblockMSL = [];
isMSL = [];
l=0;
p=0;
o=0;
for i = 1:nmouse
    idx1 = ismember(mouse_num,umouse(i));
    rule1 = rule_ordiante(idx1);
    first_sess1 = first_session_learned(idx1);
    sess1 = session_ordiante(idx1);
    correct1 = correct_trial(idx1);
    block1 = block_ordinare(idx1);
    SSL1 = SSL(idx1);
    MSL1 = MSL(idx1);
    urule = unique(rule1);
    nrule = length(urule);
    for j = 1:nrule
        idx2 = ismember(rule1,urule(j));
        first_sess2 = first_sess1(idx2);
        ii = find(first_sess2==1);
        if ~isempty(ii)
            l=l+1;
            numtr2succ(l) = ii(end);
            mnum2succ(l) = umouse(i);
        end
        sess2 = sess1(idx2);
        correct2 = correct1(idx2);
        block2 = block1(idx2);
        SSL2 = SSL1(idx2);
        MSL2 = MSL1(idx2);
        usess = unique(sess2);
        nsess = length(usess);
        for k = 1:nsess
            idx3 = ismember(sess2,usess(k));
            correct3 = correct2(idx3);
            block3 = block2(idx3);
            SSL3 = SSL2(idx3);
            MSL3 = MSL2(idx3);            
            ntr = length(correct3);
            succtr = sum(correct3);
            p_value = binopdf(succtr, ntr, prob);
            p = p+1;
            issucc(p) = p_value<0.05;
            numtr(p) = ntr;
            nmou(p) = i;
            ublock = unique(block3);
            nblock = length(ublock);
            for s = 1:nblock
                o = o+1;
                idx4 = ismember(block3,ublock(s));
                correct4 = correct3(idx4);
                SSL4 = SSL3(idx4);
                MSL4 = MSL3(idx4);                   
                succperblock(o) = sum(correct4)/length(correct4);
                issuccperblock(o) = issucc(p);
                succperblockSSL(o) = sum(correct4(SSL4))/length(correct4(SSL4));
                isSSL(o) = logical(unique(SSL4));
                succperblockMSL(o) = sum(correct4(MSL4))/length(correct4(MSL4));
                isMSL(o) = logical(unique(MSL4));                
            end
        end
    end
end

figure;
for i = 1:nmouse
    idx = ismember(mnum2succ,umouse(i));
    [ ycdf, xcdf ] = cdfcalc( numtr2succ(idx) );
    
    xx = ( xcdf * ones( 1, 2 ) )'; 
    xx = [ xx( 1 ); xx( : ) ];
    
    yy = ( ycdf( 1 : end ) * ones( 1, 2 ) )'; yy = [ yy( : ) ];
    yy( end ) = [];
    plot(xx,yy, 'Color', [0.5 0.5 0.5])
    hold on
end
[ ycdf, xcdf ] = cdfcalc( numtr2succ );

xx = ( xcdf * ones( 1, 2 ) )'; 
xx = [ xx( 1 ); xx( : ) ];

yy = ( ycdf( 1 : end ) * ones( 1, 2 ) )'; yy = [ yy( : ) ];
yy( end ) = [];
plot(xx,yy, 'Color', [0 0 0])
medi = xx( find( yy >= 0.5, 1, 'first' ) ); 
alines_local(0.5,'y','LineStyle', '--');
alines_local(medi,'x','LineStyle', '--');
set(gca, 'TickDir','out','Box','off')
axis square
ylabel('Cumulative probability')
xlabel('Testing trials to success per DC')
title('Figure S1B')

% plot figure S1C - The number of testing trials performed 
% during successful and unsuccessful sessions
figure;
for i = 1:nmouse
    idx = ismember(nmou,umouse(i)) & issucc==0;
    [ ycdf, xcdf ] = cdfcalc( numtr(idx) );
    
    xx = ( xcdf * ones( 1, 2 ) )'; 
    xx = [ xx( 1 ); xx( : ) ];
    
    yy = ( ycdf( 1 : end ) * ones( 1, 2 ) )'; yy = [ yy( : ) ];
    yy( end ) = [];
    plot(xx,yy, 'Color', [0.5 0.5 0.5])
    hold on

    idx2 = ismember(nmou,umouse(i)) & issucc==1;
    [ ycdf, xcdf ] = cdfcalc( numtr(idx2) );
    
    xx = ( xcdf * ones( 1, 2 ) )'; 
    xx = [ xx( 1 ); xx( : ) ];
    
    yy = ( ycdf( 1 : end ) * ones( 1, 2 ) )'; yy = [ yy( : ) ];
    yy( end ) = [];
    plot(xx,yy, 'Color', [1 0 0 0.25])
end

[ ycdf, xcdf ] = cdfcalc( numtr(issucc==0) );

xx = ( xcdf * ones( 1, 2 ) )'; 
xx = [ xx( 1 ); xx( : ) ];

yy = ( ycdf( 1 : end ) * ones( 1, 2 ) )'; yy = [ yy( : ) ];
yy( end ) = [];
plot(xx,yy, 'Color', [0 0 0])
medi1 = xx( find( yy >= 0.5, 1, 'first' ) ); 

[ ycdf, xcdf ] = cdfcalc( numtr(issucc==1) );

xx = ( xcdf * ones( 1, 2 ) )'; 
xx = [ xx( 1 ); xx( : ) ];

yy = ( ycdf( 1 : end ) * ones( 1, 2 ) )'; yy = [ yy( : ) ];
yy( end ) = [];
plot(xx,yy, 'Color', [1 0 0])
medi2 = xx( find( yy >= 0.5, 1, 'first' ) ); 

alines_local(0.5,'y','LineStyle', '--');
alines_local(medi1,'x','LineStyle', '--', 'Color',[0 0 0]);
alines_local(medi2,'x','LineStyle', '--', 'Color',[1 0 0]);

set(gca, 'TickDir','out','Box','off')
axis square
ylabel('Cumulative probability')
xlabel('Testing trials per session')
title('Figure S1C')
pv = utest_local(numtr(issucc==0),numtr(issucc==1))

% plot figure S1D - Success rates are higher
% during successful compared with unsuccessful sessions.
figure;
boxplot(succperblock,issuccperblock,'notch', 'on')
set(gca, 'TickDir','out','Box','off')
axis square
ylim([0 1.1])
ylabel('Success rate')
alines_local(0.5,'y','LineStyle', '--');
set(gca,'XTickLabel',{'Unsuccesful sessions';'Succesful sessions' })
title('Figure S1D')
pv = utest_local(succperblock(issuccperblock==0),succperblock(issuccperblock==1))

% plot figure S1E - Success rates are higher during SSL compared
% with MSL DCs.
figure;
nanidxMSL = isnan(succperblockMSL);
nanidxSSL = isnan(succperblockSSL);
boxplot([succperblockMSL(~nanidxMSL)'; succperblockSSL(~nanidxSSL)'],[1*isMSL(~nanidxMSL)'; 2*isSSL(~nanidxSSL)']...
    ,'notch', 'on')
set(gca, 'TickDir','out','Box','off')
axis square
ylim([0 1.1])
ylabel('Success rate')
alines_local(0.5,'y','LineStyle', '--');
set(gca,'XTickLabel',{'MSL sessions';'SSL sessions' })
title('Figure S1E')
pv = utest_local(succperblockMSL,succperblockSSL)

return  %LtL_make_FigureS1



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

