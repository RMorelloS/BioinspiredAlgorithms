%--------------------------------------------------------------
% The Firefly Algorithm
% Input parameters considered:
% Static parameters:
% The image to be segmented: im1;
% Number of thresholds: thresholds;
% Number of generations: iter;
% q value for Tsallis entropy: q
% Dinamic parameters:
% A parameter struct containing:
% Population size (pop_size) = 50;
% Upper Bound for image thresholding (UB) = 253;
% Lower Bound for image thresholding (LB) = 2.
%--------------------------------------------------------------
% Output parameters considered:
% Optimized thresholds: bests;
% Entropy value for the optimized thresholds: entropy;
% Number of generations until convergence: num_generations;
% Entropy value for each generation: generation_entropy.
%--------------------------------------------------------------
function [bests, entropy, num_generations, generation_entropy] = FF(im1,thresholds,iter,q, parameters)
num_generations = iter;
generation_entropy = zeros(iter,1);
histograma = psrGrayHistogram(im1);
NFireflies = parameters.pop_size;
% n=number of fireflies
% MaxGeneration=number of pseudo time steps
if nargin<4,   NFireflies=15; iter=30; end
n=NFireflies;  MaxGeneration=iter;

% Show info
% help firefly_simple.m
% Reset the random generator
range = [parameters.LB parameters.UB];

% ------------------------------------------------
alpha= 0.01; %50;      % Randomness 0--1 (highly random)
gamma= 1.0; %0.7;      % Absorption coefficient
delta= 0.97;      % Randomness reduction (similar to 
                % an annealing schedule)
% ------------------------------------------------

% ------------------------------------------------
% generating the initial locations of n fireflies
[fireflies, Lightn]=init_ffa(n,thresholds,range);
% Display the paths of fireflies in a figure with
% contours of the function to be optimized

% Iterations or pseudo time marching
zn=zeros(n,1);
distancias = zeros(n,n);
for i=1:MaxGeneration     %%%%% start iterations
       
% Evaluate new solutions
if q == 1
  for k=1:n  
     zn(k)= psrAvaliacaoShannon(histograma, fireflies(k,:));        
  end   
else
  for k=1:n
     zn(k)= TsallisEvaluation(histograma, q, fireflies(k,:)); 
  end
end



% Ranking the fireflies by their light intensity
[Lightn,Index]= sort(-zn);
Lightn = -Lightn;
fireflies = fireflies(Index,:);



for j=1:n
    for k=1:n
    %% outra modificacao no programa do Guilherme
    distancias(j,k) = dist(fireflies(j,:),fireflies(k,:));
    end
end

% Move all fireflies to the better locations
[fireflies]=ffa_move(fireflies,Lightn,distancias,alpha,gamma,range);
    
% Reduce randomness as iterations proceed
alpha=newalpha(alpha,delta);
melhor = fireflies(1,:);

generation_entropy(i, 1) = Lightn(1);
    
if i > 30  && std(generation_entropy(i-30:i))  < 0.01
    num_generations = i;
    break;
end

%figure(1);
%subplot(121);
%hold off;
%plot(histograma);
%title('Segmentacao');
%hold on;
%for l=1:thresholds
%line([melhor(l) melhor(l)],[min(histograma) max(histograma)],'color','red');
%end
%drawnow;
%subplot(122);
%hold off;
%plot(Lightn);
%title('Avaliação em ordem crescente');
%drawnow;

end   %%%%% end of iterations
bests=melhor;
entropy=Lightn(1);



end


% ----- All subfunctions are listed here ---------
% The initial locations of n fireflies
function [fireflies,Lightn]=init_ffa(n,thresholds,range)
  for cont=1:n
    fireflies(cont,:) = sort(randperm(254, thresholds) + 1);
  end
  Lightn=zeros(n,1);
end


% Move all fireflies toward brighter ones
function elementos=ffa_move(elementos,Lightn,distancias,alpha,gamma,range)
ni=size(elementos,1);

for i=1:ni,
% The attractiveness parameter beta=exp(-gamma*r)
    for j=1:ni
        %if Lightn(j)<=Lightn(i), % Brighter and more attractive
        if Lightn(j)>Lightn(i), % Brighter and more attractive
            beta0=1.0;     
            beta=beta0*exp(-gamma*distancias(i,j).^2);            
            vet = alpha*(randi(255,1,size(elementos,2))+1);
            if(i~=ni)
                for k=1:size(elementos,2) 
                    elementos(i,k)=(1-beta)*elementos(i,k)+beta*(elementos(j,k))+vet(1,k);   
                    elementos(i,k) = elementos(i,k)./(1+alpha);
                end
            end
        end
    end % end for j
end % end for i
end


% Reduce the randomness during iterations
function alpha=newalpha(alpha,delta)
alpha=alpha*delta;
end
%  ============== end =====================================


function ret=dist(elem1,elem2)
   ret= sqrt(sum((elem1-elem2).^2));
end

