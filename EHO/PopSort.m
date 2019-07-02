function [Population, indices] = PopSort(Population, range)
% Sort the population members from best to worst
popsize = length(Population);
Cost = zeros(1, popsize);
indices = zeros(1, popsize);
for i = 1 : popsize
    Cost(i) = Population(i).cost;
end
[Cost, indices] = sort(Cost, 2, 'descend');

Chroms = zeros(popsize, length(Population(1).chrom));
for i = 1 : popsize
    vetor = sort(Population(indices(i)).chrom);
    Chroms(i, :) = vetor;
end
for i = 1 : popsize
    Population(i).chrom = Chroms(i, :);
    Population(i).cost = Cost(i);
end


