
%% Elephant Herding Optimization (EHO)
% Author: Gai-Ge Wang
% Email: gaigewang@163.com
%             gaigewang@gmail.com

% Main paper:
% Gai-Ge Wang, Suash Deb, and Leandro dos Santos Coelho, Elephant Herding Optimization.
% In: 2015 3rd International Symposium on Computational
% and Business Intelligence (ISCBI 2015), 
% Bali, Indonesia, December 7-8, 2015. IEEE, pp ****


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%% Notes:
% Different run may generate different solutions, this is determined by
% the the nature of metaheuristic algorithms.
%%

function [MinCost] = EHO(ProblemFunction, DisplayFlag, RandSeed)

%  Elephant Herding Optimization (EHO) software for minimizing a general function
% The fixed generation is considered as termination condition.

% INPUTS: ProblemFunction is the handle of the function that returns
%         the handles of the initialization, cost, and feasibility functions.
%         DisplayFlag = true or false, whether or not to display and plot results.
%         ProbFlag = true or false, whether or not to use probabilities to update emigration rates.
%         RandSeed = random number seed
% OUTPUTS: MinCost = array of best solution, one element for each generation
%         
% CAVEAT: The "ClearDups" function that is called below replaces duplicates with randomly-generated
%         individuals, but it does not then recalculate the cost of the replaced individuals.

tic

if ~exist('ProblemFunction', 'var')
    ProblemFunction = @Ackley;
end
if ~exist('DisplayFlag', 'var')
    DisplayFlag = true;
end
if ~exist('RandSeed', 'var')
    RandSeed = round(sum(100*clock));
end

[OPTIONS, MinCost, AvgCost, InitFunction, CostFunction, FeasibleFunction, ...
    MaxParValue, MinParValue, Population] = Init(DisplayFlag, ProblemFunction, RandSeed);

% % % % % % % % % % % %             Initial parameter setting          % % % % % % % % % % % %%%%
%% Initial parameter setting
Keep = 2; % elitism parameter: how many of the best elephants to keep from one generation to the next
numElephantInEachClan = OPTIONS.numElephantInEachClan*ones(1, OPTIONS.numClan);
dim = OPTIONS.numVar;
alpha = 0.5;
beta = 0.1;
% % % % % % % % % % % %       End of Initial parameter setting       % % % % % % % % % % % %%
%%

% % % % % % % % % % % %             Begin the optimization loop        % % % % % % % % % %%%%
% Begin the optimization loop
for GenIndex = 1 : OPTIONS.Maxgen
    
    % % % % % % % % % % % %            Elitism Strategy           % % % % % % % % % % % %%%%%
    %% Save the best elephants in a temporary array.
    for j = 1 : Keep
        chromKeep(j,:) = Population(j).chrom;
        costKeep(j) = Population(j).cost;
    end  %end for j
    % % % % % % % % % % % %       End of  Elitism Strategy      % % % % % % % % % % % %%%%
    %%
    
    % % % % % % % % % % % %    Divide the whole population into some clans % % % %%%% % % %%%
    %% Divide the whole elephant population into some clans
    % according to their fitness.
    %% 分解成若干个clan
    j =1;
    popindex = 1;
    while popindex <= OPTIONS.popsize
        
        for cindex = 1 : OPTIONS.numClan
            Clan{cindex}(j) = Population(popindex);
            popindex = popindex + 1;
        end % end for cindex
        
        j = j+1;
    end  % end for popindex
    % % % % % % % % % % %    End of Divide the whole population into some clans  % % %%% % % % % %
    %%
    
    % % % % % % % % % % % %%            Clan updating operator          % % % % % % % % % % % %%%%
    %% Clan updating operator
    j = 1;
    popindex =1;
    while popindex <= OPTIONS.popsize
        
        for cindex = 1 : OPTIONS.numClan
            ClanCenter =  CaculateClanCenter(OPTIONS, Clan, cindex);
            NewClan{cindex}(j).chrom = Clan{cindex}(j).chrom ...
                + alpha*(Clan{cindex}(1).chrom - Clan{cindex}(j).chrom).* rand(1,dim);
            
            if sum( NewClan{cindex}(j).chrom - Clan{cindex}(j).chrom) == 0
                NewClan{cindex}(j).chrom = beta * ClanCenter ;
            end
            
            popindex = popindex + 1;
        end % end for cindex
        
        j = j+1;
    end % end for popindex
    % % % % % % % % % % % %%%       End of Clan updating operator      % % % % % % % % % % %%%
    %%
    
    % % % % % % % % % % % %%            Separating operator          % % % % % % % % % % % %%%%
    %% Separating operator
    for cindex = 1 : OPTIONS.numClan
        NewClan{cindex}(end).chrom = MinParValue...
            + (MaxParValue - MinParValue + 1) .* rand(1,OPTIONS.numVar);
    end % end for cindex
    % % % % % % % % % % % %%%       End of Separating operator      % % % % % % % % % % % %%%
    %%
    
    % % % % % % % % % % %             Evaluate NewClan          % % % % % % % % % % % %% % % %
    %% Evaluate NewClan
    SavePopSize = OPTIONS.popsize;
    for i=1:OPTIONS.numClan
        OPTIONS.popsize = numElephantInEachClan(i);
        % Make sure the population does not have duplicates.
        NewClan{i} = ClearDups(NewClan{i}, MaxParValue, MinParValue);
        % Make sure each individual is legal.
        NewClan{i} = FeasibleFunction(OPTIONS, NewClan{i});
        % Calculate cost
        NewClan{i} = CostFunction(OPTIONS, NewClan{i});
        % Sort from best to worst
        NewClan{i} = PopSort(NewClan{i});
    end % end for i
    OPTIONS.popsize = SavePopSize;
    % % % % % % % % % % % %       End of Evaluate NewClan     % % % % % % % % % % % %%
    %%
    
    % % % % % % %  Combine two subpopulations into one and rank monarch butterflis       % % % % % %
    %% Combine Population1 with Population2 to generate a new Population
    Population = CombineClan(OPTIONS, NewClan);
    % Sort from best to worst
    Population = PopSort(Population);
    % % % % % %     End of Combine two subpopulations into one and rank monarch butterflis  % %% % %
    %%
    
    % % % % % % % % % % % %            Elitism Strategy          % % % % % % % % % % % %%% %% %
    %% Replace the worst with the previous generation's elites.
    n = length(Population);
    for k3 = 1 : Keep
        Population(n-k3+1).chrom = chromKeep(k3,:);
        Population(n-k3+1).cost = costKeep(k3);
    end % end for k3
    % % % % % % % % % % % %     End of  Elitism Strategy      % % % % % % % % % % % %%% %% %
    %%
    
    % % % % % % % % % %           Precess and output the results          % % % % % % % % % % % %%%
    % Sort from best to worst
    Population = PopSort(Population);
    % Compute the average cost
    [AverageCost, nLegal] = ComputeAveCost(Population);
    % Display info to screen
    MinCost = [MinCost Population(1).cost];
    AvgCost = [AvgCost AverageCost];
    if DisplayFlag
        disp(['The best and mean of Generation # ', num2str(GenIndex), ' are ',...
            num2str(MinCost(end)), ' and ', num2str(AvgCost(end))]);
    end
    % % % % % % % % % % %    End of Precess and output the results     %%%%%%%%%% %% %
    %%
    
    
end % end for GenIndex
Conclude1(DisplayFlag, OPTIONS, Population, nLegal, MinCost, AvgCost);

save Population;

toc

% % % % % % % % % %     End of Elephant Herding Optimization implementation     %%%% %% %
%%
