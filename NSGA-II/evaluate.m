function [pop, state] = evaluate(opt, pop, state, varargin)
% Function: [pop, state] = evaluate(opt, pop, state, varargin)
% Description: Evaluate the objective functions of each individual in the
%   population.
%
%         LSSSSWC, NWPU
%    Revision: 1.0  Data: 2011-04-20
%*************************************************************************

N = length(pop);
allTime = zeros(N, 1);  % allTime : use to calculate average evaluation times

%*************************************************************************
% Evaluate objective function in parallel
%*************************************************************************
if( strcmpi(opt.useParallel, 'yes') == 1 )
%     curPoolInfo = gcp('nocreate');
%     
%     % Check if there is a parpool open or not
%     [Pooldatasize, ~]=size(curPoolInfo);
%     if Pooldatasize == 0
%         curPoolsize=0;
%     else
%         curPoolsize = curPoolInfo.NumWorkers;
%     end
% 
%     % There isn't opened worker process
%     if(curPoolsize == 0)
%         if(opt.poolsize == 0)
%             parpool open local;
%         else
%             parpool(opt.poolsize);
%         end
%     % Close and recreate worker process
%     else
%         if(opt.poolsize ~= curPoolsize)
%             delete(gcp);
%             parpool(opt.poolsize);
%         end
%     end

    parfor i = 1:N
        fprintf('\nEvaluating the objective function... Generation: %d / %d, Individual: %d / %d \n',state.currentGen, opt.maxGen, i, N);
        [pop(i), allTime(i)] = evalIndividual(pop(i), opt.objfun, varargin{:});
    end

%*************************************************************************
% Evaluate objective function in serial
%*************************************************************************
else
    for i = 1:N
        fprintf('\nEvaluating the objective function... Generation: %d / %d, Individual: %d / %d \n',state.currentGen, opt.maxGen, i, N);
        [pop(i), allTime(i)] = evalIndividual(pop(i), opt.objfun, varargin{:});
    end
end

%*************************************************************************
% Statistics
%*************************************************************************
state.avgEvalTime   = sum(allTime) / length(allTime);
state.evaluateCount = state.evaluateCount + length(pop);




function [indi, evalTime] = evalIndividual(indi, objfun, varargin)
% Function: [indi, evalTime] = evalIndividual(indi, objfun, varargin)
% Description: Evaluate one objective function.
%
%         LSSSSWC, NWPU
%    Revision: 1.1  Data: 2011-07-25
%*************************************************************************

tStart = tic;
[y, cons] = objfun( indi.var, varargin{:} );
evalTime = toc(tStart);

% Save the objective values and constraint violations
indi.obj = y;

if( ~isempty(indi.cons) )
    idx = find( cons );
    indi.cons=cons;
    if( ~isempty(idx) )
        indi.nViol = length(idx);
        indi.violSum = sum( abs(cons) );
    else
        indi.nViol = 0;
        indi.violSum = 0;
    end

end


