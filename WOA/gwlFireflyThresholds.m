
%
% histograma: histograma em tons de cinza de uma imagem 
% thresholds: 1, 2, 3, 4, 5 ...
% NFireflies: 50
% iter: 100
% METHOD: 'TE'
% bests: é um vetor de limiares
%
%
function bests = gwlFireflyThresholds(histograma,thresholds,NFireflies,iter,q,METHOD)

% n=number of fireflies
% MaxGeneration=number of pseudo time steps
if nargin<4,   NFireflies=15; iter=30; end
n=NFireflies;  MaxGeneration=iter;

% Show info
% help firefly_simple.m
%rand('state',0);  % Reset the random generator
range = [2 253];

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
for i=1:MaxGeneration,     %%%%% start iterations
       
% Evaluate new solutions
switch METHOD
    case {'CE'}
      for k=1:n  
         zn(k)= psrAvaliacaoCrossEntropy(histograma, fireflies(k,:),1);      
      end
    case {'CEM'}
      for k=1:n  
         zn(k)= psrAvaliacaoCrossEntropy(histograma, fireflies(k,:),6);      
      end
    case {'SE'}
      for k=1:n  
         zn(k)= psrAvaliacaoShannon(histograma, fireflies(k,:));        
      end   
    case {'TE'}
      for k=1:n  
         zn(k)= psrAvaliacaoTsallis(histograma, q, fireflies(k,:)); 
      end
    otherwise
    disp('error in METHOD of evaluation');
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

end


% ----- All subfunctions are listed here ---------
% The initial locations of n fireflies
function [fireflies,Lightn]=init_ffa(n,thresholds,range)

  fireflies = sort(randi(253,n, thresholds)+1,2);
  Lightn=zeros(n,1);
  fireflies=validaElementos(fireflies,range);

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
elementos=validaElementos(elementos,range);

end


% Reduce the randomness during iterations
function alpha=newalpha(alpha,delta)
alpha=alpha*delta;
end

% Make sure the fireflies are within the range
function elementos=validaElementos(elementos, range)
elementos=round(elementos);
elementos=sort(elementos,2);

    for i=1:size(elementos,1)
        
        %evita que os thresholds tenham o mesmo valor
        if(size(elementos(i,:),2)~=size(unique(elementos(i,:)),2))
            elementos(i,:) = sort(randi(253,1,size(elementos,2))+1);
        end
       
        for c=1:size(elementos,2)
            
            if(elementos(i,c)>range(2))%%% pra que isso?
                %elementos(i,c) = randi(252)+1;
                elementos(i,:)=elementos(i,:)/elementos(i,c)*255;
                elementos(i,:)=round(elementos(i,:));
            elseif elementos(c)<range(1)
                elementos(i,c) = randi(252)+1; 
            end
            
        end
    end
end
%  ============== end =====================================


function ret=dist(elem1,elem2)
   ret= sqrt(sum((elem1-elem2).^2));
end

