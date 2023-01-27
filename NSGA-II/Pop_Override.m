function pop= Pop_Override(opt, pop, varargin)
% Function: pop = Pop_Override(opt, pop, varargin)
% Description: Load population from exist variables. Does not reject based
% on number of constraints from previous results

    popsize_new=opt.popsize;

    variable=opt.initfun{2};
    
    popsize_old=size(variable, 1);

    numb_var_sim=opt.numVar;

    numb_var_supp=size(variable, 2);

    if numb_var_sim ~= numb_var_supp
        error('NSGA2:OptModelError', 'Number of variables supplied does not equal what is requested');
    end

    pop_size=opt.popsize;

    variable=opt.initfun{2};
    
    if popsize_old >= popsize_new
        for j = 1:pop_size
            pop(j).var = variable(j, :);
        end
    else
        for j = 1:popsize_old
            pop(j).var = variable(j, :);
        end
        pop(popsize_old+1:end) = initpopUniform(opt, pop(popsize_old+1:end));
    end
    
    
    
    function pop = initpopUniform(opt, pop)
    % Function: pop = initpopUniform(opt, pop)
    % Description: Initialize population using random number
    %
    %    Copyright 2011 by LSSSSWC
    %    Revision: 1.0  Data: 2011-07-01
    %*************************************************************************

        nVar = opt.numVar;
        type = opt.vartype;

        lb = opt.lb;
        ub = opt.ub;

        popsize = length(pop);
        for i = 1:popsize
            var = lb + rand(1, nVar) .* (ub-lb);

            % if desing variable is integer, round to the nearest integer
            for v = 1:nVar
                if( type(v) == 2)
                    var(v) = round(var(v));
                end
            end

            % limit in the lower and upper bound
            var = varlimit(var, lb, ub);

            pop(i).var = var;

        end
    end

end