
function [bests, cost, numGeracoes, melhores_avaliacoes] = EHO(image, thresholds, geracoes, q, parametrosEHO)

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

[~, Population,  numGeracoes, melhores_avaliacoes] = EHO_FEs_V2(ProblemFunction, DisplayFlag, RandSeed, image, geracoes, q, thresholds, parametrosEHO);
for i=1:100     %matriz de indivíduos
  for j=1:nl
    R(i,j) = Population(i).chrom(:,j);
  end
  
  R(i,nl+1) = Population(i).cost;
end
cost = R(1,nl+1);
bests = round(R(1,1:nl));

end





