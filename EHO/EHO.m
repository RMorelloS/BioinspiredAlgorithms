%--------------------------------------------------------------
% The Elephant Herding Optimization Algorithm
% Input parameters considered:
% Static parameters:
% The image to be segmented: im;
% Number of thresholds: thresholds;
% Number of generations: Max_iter;
% q value for Tsallis entropy: q
% Dinamic parameters:
% A parameter struct containing:
% Population size (pop_size) = 100;
% Number of clans (numClan) = 5;
% Elitism (Keep) = 2;
% Alpha (alpha) = 0.5;
% Beta (beta) = 0.1;
% Upper Bound for image thresholding (UB) = 253;
% Lower Bound for image thresholding (LB) = 2.
%--------------------------------------------------------------
% Output parameters considered:
% Optimized thresholds: best;
% Entropy value for the optimized thresholds: cost;
% Number of generations until convergence: num_generations;
% Entropy value for each generation: generation_entropy.
%--------------------------------------------------------------
function [bests, cost, num_generations, generation_entropy] = EHO(image, thresholds, Max_iter, q, parameters)

if ~exist('ProblemFunction', 'var')
    ProblemFunction = @Ackley;
end
if ~exist('DisplayFlag', 'var')
    DisplayFlag = true;
end
if ~exist('RandSeed', 'var')
    RandSeed = round(sum(100*clock));
end

nl = thresholds; % number of variables in each population member;

[~, Population,  num_generations, generation_entropy] = EHO_FEs_V2(ProblemFunction, DisplayFlag, RandSeed, image, Max_iter, q, thresholds, parameters);
for i=1:parameters.pop_size     %matriz de indivíduos
  for j=1:nl
    R(i,j) = Population(i).chrom(:,j);
  end
  
  R(i,nl+1) = Population(i).cost;
end
cost = R(1,nl+1);
bests = round(R(1,1:nl));

end





