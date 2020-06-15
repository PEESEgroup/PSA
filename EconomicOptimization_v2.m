clc
format long
parpool('local', 12);

addpath('CycleSteps')
addpath('GA_files')

load('Params')

N = 30 ;
type = 'EconomicEvaluation' ;

for i = 12:12
    
% load parameters
IsothermParams     = IsothermPar(i, :) ;    
material_propertry = SimParam(i, :)    ;

material    = {}                 ;
material{1} = material_propertry ;
material{2} = IsothermParams     ;

Function = @(x) PSACycleSimulation( x, material, type, N ) ; % Function to simulate the PSA cycle

% initial variables
[~, vars] = sortt(loadpopfile('Ec.txt'));

options            = nsgaopt();                          % create default options structure
options.popsize    = 60;                                 % populaion size
options.outputfile = 'Ec_2.txt';
options.maxGen     = 120;                                % max generation

options.vartype    = [1, 1, 1, 1, 1, 1] ;

options.initfun={@Pop_Override, vars}   ;                % Supply variables from previous results

options.numObj  = 2 ;                                    % number of objectives
options.numVar  = 6 ;                                    % number of design variables
options.numCons = 3 ;                                    % number of constraints
options.lb = [1e5,  10, 0.01, 0.1, 0, 1e4]   ;           % lower bound of x
options.ub = [10e5, 1000, 0.99, 2, 1, 5e4]   ;           % upper bound of x
options.nameObj = {'-productivity','energy'} ;           % the objective names are showed in GUI window.
options.objfun  = Function                   ;           % objective function handle

options.useParallel = 'yes' ;                            % parallel computation is non-essential here
options.poolsize     = 12   ;                            % number of worker processes

result2 = nsga2(options)    ;                            % begin the optimization!

end

