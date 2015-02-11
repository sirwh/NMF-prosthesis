%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: A slightly customed NMF function based on the original matlab NMF
% algorithm "nnmf.m". 
%--------------------------------------------------------------------------------------------------------------------------------
% Used in: Custom_NMF.m
%--------------------------------------------------------------------------------------------------------------------------------
% Update : 2015/2/12 by Jiaxin MA @Kyoto Univ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [wbest,hbest,normbest] = custom_nnmf(a,k,varargin)

% Check required arguments
error(nargchk(2,Inf,nargin,'struct'))
[n,m] = size(a);
if ~isscalar(k) || ~isnumeric(k) || k<1 || k>min(m,n) || k~=round(k)
    error(message('stats:nnmf:BadK'));
end

% Process optional arguments
pnames = {'algorithm' 'w0' 'h0' 'replicates' 'options'};
dflts  = {'als'       []   []   1            []       };
[eid,emsg,alg,w0,h0,tries,options] = ...
        internal.stats.getargs(pnames,dflts,varargin{:});
if ~isempty(eid)
    error(sprintf('stats:nnmf:%s',eid),emsg)
end
    
% Check optional arguments
% alg = statgetkeyword(alg,{'mult' 'als'},false,'ALGORITHM','stats:nnmf:BadAlg');
ismult = strncmp('mult',alg,numel(alg));
checkmatrices(a,w0,h0,k);
if ~isscalar(tries) || ~isnumeric(tries) || tries<1 || tries~=round(tries)
    error(message('stats:nnmf:BadReplicates'));
end

defaultopt = statset('nnmf');
tolx = statget(options,'TolX',defaultopt,'fast');
tolfun = statget(options,'TolFun',defaultopt,'fast');
maxiter = statget(options,'MaxIter',defaultopt,'fast');
dispopt = statget(options,'Display',defaultopt,'fast');

dispnum = strmatch(lower(dispopt), {'off','notify','final','iter'}) - 1;

[useParallel, RNGscheme, poolsz] = ...
    internal.stats.parallel.processParallelAndStreamOptions(options,true);

usePool = useParallel && poolsz>0;

% Modified part
unchange = [];
if ~isempty(w0) && isempty(h0)
    unchange = 'w0';
elseif ~isempty(h0) && isempty(w0)
    unchange = 'h0';
end

% Special case, if K is full rank we know the answer
if isempty(w0) && isempty(h0)
    if k==m
        w0 = a;
        h0 = eye(k);
    elseif k==n
        w0 = eye(k);
        h0 = a;
    end
end


% Define the function that will perform one iteration of the
% loop inside smartFor
loopbody = @loopBody;

% Suppress undesired warnings.
if usePool
    % On workers and client
    pctRunOnAll internal.stats.parallel.muteParallelStore('rankDeficientMatrix', ...
        warning('off','MATLAB:rankDeficientMatrix') );
else
    % On client
    ws = warning('off','MATLAB:rankDeficientMatrix');
end

% Prepare for in-progress 
if dispnum > 1 % 'iter' or 'final'
    if usePool
       % If we are running on a matlabpool, each worker will generate
       % a separate periodic report.  Before starting the loop, we
       % seed the matlabpool so that each worker will have an
       % identifying label (eg, index) for its report.
       internal.stats.parallel.distributeToPool( ...
           'workerID', num2cell(1:poolsz) );

       % Periodic reports behave differently in parallel than they do
       % in serial computation (which is the baseline).
       % We advise the user of the difference.
       warning('stats:nnmf:displayParallel',['When using the display option in'...
         ,' ','parallel, the "replicate" indices will be out of order.'...
         ,' ','The "worker" column shows which worker performed each replicate.']) 
       fprintf('    worker\t      rep\t   iteration\t     rms resid\t     |delta x|\n' );
    else
       if useParallel
          warning(message('stats:nnmf:displayParallel'));
       end
       fprintf('    rep\t   iteration\t   rms resid\t  |delta x|\n');
    end
end

try
    whbest = internal.stats.parallel.smartForReduce(...
        tries, loopbody, useParallel, RNGscheme, 'argmin');
catch ME
    % Revert warning setting for rankDeficientMatrix to value prior to nnmf.
    if usePool
        % On workers and on client
        pctRunOnAll warning(internal.stats.parallel.statParallelStore('rankDeficientMatrix').state,'MATLAB:rankDeficientMatrix');
    else
        % On client
        warning(ws);
    end
    rethrow(ME);
end

normbest = whbest{1};
wbest = whbest{3};
hbest = whbest{4};
% whbest{2} contains the iteration chosen for the best factorization,
% but it has no meaning except as a "reproducible" tie-breaker, and
% is not supplied as a return value.
    
if dispnum > 1   % 'final' or 'iter'
    fprintf('Final root mean square residual = %g\n',normbest);
end

% Revert warning setting for rankDeficientMatrix to value prior to nnmf.
if usePool
    % On workers and on client
    pctRunOnAll warning(internal.stats.parallel.statParallelStore('rankDeficientMatrix').state,'MATLAB:rankDeficientMatrix');
else
    % On client
    warning(ws);
end

if normbest==Inf
    error(message('stats:nnmf:NoSolution'))
end

% Put the outputs in a standard form - first normalize h
hlen = sqrt(sum(hbest.^2,2));
if any(hlen==0)
    warning(message('stats:nnmf:LowRank', k - sum( hlen==0 ), k));
    hlen(hlen==0) = 1;
end
% wbest = bsxfun(@times,wbest,hlen');
% hbest = bsxfun(@times,hbest,1./hlen);

% Then order by w
[~,idx] = sort(sum(wbest.^2,1),'descend');
% wbest = wbest(:,idx);
% hbest = hbest(idx,:);

% ---- Nested functions ----

    function cellout = loopBody(iter,S)
        if isempty(S)
            S = RandStream.getDefaultStream;
        end

        % whtry is a "temporary variable" and hence needs to be
        % reinitialized at start of each loop.  
        whtry = cell(4,1); % whtry{1} = norm of error
                           % whtry{3} = w
                           % whtry{4} = h

        % Get random starting values if required
        if( ~isempty(w0) && iter ==1 ) 
            whtry{3} = w0;
        else
            whtry{3} = rand(S,n,k); 
        end
        if( ~isempty(h0) && iter ==1 )
            whtry{4} = h0;
        else
            whtry{4} = rand(S,k,m); 
        end

        % Perform a factorization
        [whtry{3},whtry{4},whtry{1}] = ...
            nnmf1(a,whtry{3},whtry{4},ismult,maxiter,tolfun,tolx,...
                                  dispnum,iter,usePool,unchange);
        whtry{2} = iter;

        cellout = whtry;
    end
end  % of nnmf

% -------------------
function [w,h,dnorm] = nnmf1(a,w0,h0,ismult,maxiter,tolfun,tolx,...
                                   dispnum,repnum,usePool,unchange)
% Single non-negative matrix factorization
nm = numel(a);
sqrteps = sqrt(eps);


% Display progress.  For parallel computing, the replicate number will be
% displayed under the worker performing the replicate.
if dispnum>1 % 'final' or 'iter' 
    if usePool 
        labindx = internal.stats.parallel.workerGetValue('workerID');
        dispfmt = '%8d\t%8d\t%8d\t%14g\t%14g\n';
    else
        dispfmt = '%7d\t%8d\t%12g\t%12g\n';
    end   
end    

% Modified part
for j=1:maxiter
    if ismult
        % Multiplicative update formula
        numer = w0'*a;
        h = h0 .* (numer ./ ((w0'*w0)*h0 + eps(numer)));
        if strcmp(unchange,'h0')
            h=h0;
        end
        numer = a*h';
        w = w0 .* (numer ./ (w0*(h*h') + eps(numer)));
        if strcmp(unchange,'w0')
            w=w0;
        end
    else
        % Alternating least squares
        h = max(0, w0\a);
        w = max(0, a/h);
    end

    
    % Get norm of difference and max change in factors
    d = a - w*h;
    dnorm = sqrt(sum(sum(d.^2))/nm);
    dw = max(max(abs(w-w0) / (sqrteps+max(max(abs(w0))))));
    dh = max(max(abs(h-h0) / (sqrteps+max(max(abs(h0))))));
    delta = max(dw,dh);
%     if strcmp(unchange,'h0')
%         delta=dw;
%     elseif strcmp(unchange,'w0')
%         delta=dh;
%     end
    
    % Check for convergence
    if j>1
        if delta <= tolx
            break;
        elseif dnorm0-dnorm <= tolfun*max(1,dnorm0)
            break;
        elseif j==maxiter
            break
        end
    end

    if dispnum>2 % 'iter'
       if usePool 
           fprintf(dispfmt,labindx,repnum,j,dnorm,delta)
       else
           fprintf(dispfmt,repnum,j,dnorm,delta);
       end
    end

    % Remember previous iteration results
    dnorm0 = dnorm;
    w0 = w;
    h0 = h;
    
end

if dispnum>1   % 'final' or 'iter'
    if usePool 
       fprintf(dispfmt,labindx,repnum,j,dnorm,delta)
   else
       fprintf(dispfmt,repnum,j,dnorm,delta);
   end
end

end

% ---------------------------
function checkmatrices(a,w,h,k)
% check for non-negative matrices of the proper size

if ndims(a)~=2 || ~isnumeric(a) || ~isreal(a) || any(any(a<0)) || any(any(~isfinite(a)))
    error(message('stats:nnmf:BadA'))
end
[n,m] = size(a);
if ~isempty(w)
    if ndims(w)~=2 || ~isnumeric(w)|| ~isreal(w) || any(any(w<0)) || any(any(~isfinite(w)))
        error('stats:nnmf:BadW','W must be a matrix of non-negative values.')
    elseif ~isequal(size(w),[n k])
        error('stats:nnmf:BadW',...
              'The size of the W matrix must be %d-by-%d.',n,k)
    end
end
if ~isempty(h)
    if ndims(h)~=2 || ~isnumeric(h)|| ~isreal(h) || any(any(h<0)) || any(any(~isfinite(h)))
        error('stats:nnmf:BadH','H must be a matrix of non-negative values.')
    elseif ~isequal(size(h),[k m])
        error('stats:nnmf:BadH',...
              'The size of the H matrix must be %d-by-%d.',k,m)
    end
end
end % checkmatrices

