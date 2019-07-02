%_________________________________________________________________________%
%  Whale Optimization Algorithm (WOA) source codes demo 1.0               %
%                                                                         %
%  Developed in MATLAB R2011b(7.13)                                       %
%                                                                         %
%  Author and programmer: Seyedali Mirjalili                              %
%                                                                         %
%         e-Mail: ali.mirjalili@gmail.com                                 %
%                 seyedali.mirjalili@griffithuni.edu.au                   %
%                                                                         %
%       Homepage: http://www.alimirjalili.com                             %
%                                                                         %
%   Main paper: S. Mirjalili, A. Lewis                                    %
%               The Whale Optimization Algorithm,                         %
%               Advances in Engineering Software , in press,              %
%               DOI: http://dx.doi.org/10.1016/j.advengsoft.2016.01.008   %
%                                                                         %
%_________________________________________________________________________%


% The Whale Optimization Algorithm
function [Leader_pos, Leader_score, numGeracoes, melhores_avaliacoes]=WOA(im1, dim, Max_iter, q, parametros)
numGeracoes = Max_iter;
melhores_avaliacoes = zeros(Max_iter,1);
SearchAgents_no = parametros.nBaleias;

functionName = 'F0';
[fobj] = Get_Functions_detailsWOA(functionName);
H = psrGrayHistogram(im1);
% initialize position vector and score for the leader
Leader_pos=zeros(1,dim);
Leader_score=-inf; %change this to -inf for maximization problems
ub = parametros.UB;
lb = parametros.LB;

%Initialize the positions of search agents
Positions=initializationWOA(SearchAgents_no,dim,ub,lb);
Positions = sort(Positions')';
Convergence_curve=zeros(1,Max_iter);

t=0;% Loop counter

% Main loop

while t<Max_iter
    for i=1:size(Positions,1)
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        % Calculate objective function for each search agent
        fitness=fobj(Positions(i,:), H, q);
        % Update the leader
        if fitness>Leader_score % Change this to > for maximization problem
            Leader_score=fitness; % Update alpha
            Leader_pos=Positions(i,:);
        end
    end
    a=2-t*((2)/Max_iter); % a decreases linearly fron 2 to 0 in Eq. (2.3)
    
    % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2=-1+t*((-1)/Max_iter);
    
    % Update the Position of search agents 
    for i=1:size(Positions,1)
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]
        
        A=2*a*r1-a;  % Eq. (2.3) in the paper
        C=2*r2;      % Eq. (2.4) in the paper
        
        
        b=1;               %  parameters in Eq. (2.5)
        l=(a2-1)*rand+1;   %  parameters in Eq. (2.5)
        
        p = rand();        % p in Eq. (2.6)
        
        for j=1:size(Positions,2)
            
            if p<0.5   
                if abs(A)>=1
                    rand_leader_index = floor(SearchAgents_no*rand()+1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand=abs(C*X_rand(j)-Positions(i,j)); % Eq. (2.7)
                    Positions(i,j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8)
                    
                elseif abs(A)<1
                    D_Leader=abs(C*Leader_pos(j)-Positions(i,j)); % Eq. (2.1)
                    Positions(i,j)=Leader_pos(j)-A*D_Leader;      % Eq. (2.2)
                end
                %Positions = validaElementos(Positions, [2 253]);
            elseif p>=0.5
              
                distance2Leader=abs(Leader_pos(j)-Positions(i,j));
                % Eq. (2.5)
                Positions(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos(j);
                %Positions = validaElementos(Positions, [2 253]);
            end
            
        end

    end
    t=t+1;
    Convergence_curve(t)=Leader_score;
    [t Leader_score];
    melhores_avaliacoes(t+1, 1) = Leader_score;
    
    if t > 30  && std(melhores_avaliacoes(t-30:t))  < 0.01
       numGeracoes = t;
       break;
    end
end



