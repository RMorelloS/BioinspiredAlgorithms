function [OPTIONS, MinCost, AvgCost, InitFunction, CostFunction, FeasibleFunction, ...
    MaxParValue, MinParValue, Population] = Init(DisplayFlag, ProblemFunction, RandSeed, Histogram, geracoes, q, thresholds, parametrosEHO, range)

% Initialize population-based optimization software.
% WARNING: some of the optimization routines will not work if population size is odd.
OPTIONS.popsize = parametrosEHO.popsize; %50; % total population size
OPTIONS.Maxgen = geracoes; %200; %50; % generation count limit
OPTIONS.numVar = thresholds; % original 20; % number of vriables in each population member
OPTIONS.numClan = parametrosEHO.numClan; % number of Function Evaluations (FEs)
OPTIONS.numElephantInEachClan = OPTIONS.popsize/OPTIONS.numClan;
OPTIONS.MaxFEs = 0.3E4; % number of Function Evaluations (FEs)

if ~exist('RandSeed', 'var')
    RandSeed = round(sum(100*clock));
end
rand('state', RandSeed); % initialize random number generator
if DisplayFlag
  %  disp(['random # seed = ', num2str(RandSeed)]);
end

% Get the addresses of the initialization, cost, and feasibility functions.
[InitFunction, CostFunction, FeasibleFunction] = ProblemFunction();
% Initialize the population.
[MaxParValue, MinParValue, Population, OPTIONS] = InitFunction(OPTIONS);
% Make sure the population does not have duplicates. 
Population = ClearDups(Population, MaxParValue, MinParValue);
% Compute cost of each individual  
Population = CostFunction(OPTIONS, Population, Histogram, q);
% Sort the population from most fit to least fit
Population = PopSort(Population, range);
% Compute the average cost
AverageCost = ComputeAveCost(Population);
% Display info to screen
MinCost = [Population(1).cost];
AvgCost = [AverageCost];
%if DisplayFlag 
   % disp(['The best and mean of Generation # 0 are ', num2str(MinCost(end)), ' and ', num2str(AvgCost(end))]);
%end

return;
