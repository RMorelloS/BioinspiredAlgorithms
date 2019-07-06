% Krill Herd Algorithm V 1.1
% By Amir H. Gandomi (a.h.gandomi@gmail.com)

% Main paper:
% Gandomi A.H., Alavi A.H., Krill Herd: A New Bio-Inspired Optimization Algorithm.
% Communications in Nonlinear Science and Numerical Simulation, 
% 2012, 17 (12):4831-4845.
% DOI: 10.1016/j.cnsns.2012.05.010

% Boundary Constraint Handling Scheme used in this code:
% Gandomi A.H., Yang X.S., "Evolutionary Boundary Constraint Handling Scheme"
% Neural Computing & Applications, 2012, 21(6):1449-1462.
% DOI: 10.1007/s00521-012-1069-0
%--------------------------------------------------------------
% The Krill Herd Algorithm
% Evaluation function changed to Tsallis entropy. Changed by 
% Ricardo Morello Santos, Prof. Guilherme Wachs Lopes, 
% Prof. Nilson Saito and Prof. Paulo Sérgio Rodrigues
%--------------------------------------------------------------
% Input parameters considered:
% Static parameters:
% The image to be segmented: im;
% Number of thresholds: thresholds;
% Number of generations: Max_iter;
% q value for Tsallis entropy: q
% Dinamic parameters:
% A parameter struct containing:
% Population size (pop_size) = 40;
% Upper Bound for image thresholding (UB) = 253;
% Lower Bound for image thresholding (LB) = 2.
%--------------------------------------------------------------
% Output parameters considered:
% Optimized thresholds: Lim;
% Entropy value for the optimized thresholds: entropy;
% Number of generations until convergence: num_generations;
% Entropy value for each generation: generation_entropy.
%--------------------------------------------------------------


function [Lim, entropy, num_generations, generation_entropy] = KH(im, thresholds, Max_iter, q, parameters)
num_generations = Max_iter;
generation_entropy = zeros(Max_iter,1);
%% Initial Parameter Setting
NR = 1;  % original 10                    % Number if Runs
pop_size = parameters.pop_size; %original = 25;  % Number if Krills
MI = Max_iter; %original 200   % Maximum Iteration
C_flag = 1;                               % Crossover flag [Yes=1]

% Bounds (Normalize search space in case of highly imbalanced search space)
dim = thresholds; % original 10 % numero de limiares
UB = 253*ones(1,dim);
LB = 2*ones(1,dim);
entropy = 0;

NP = length(LB); % Number if Parameter(s)
Dt = mean(abs(UB-LB))/2; % Scale Factor

F = zeros(NP,pop_size);D = zeros(1,pop_size);N = zeros(NP,pop_size); %R = zeros(NP,pop_size);
Vf = 0.02; % original 0.02; 
Dmax = 0.005; Nmax = 0.01; Sr = 0;
xmin = 0; xmax = 0.08; ymin = 0; ymax = 0.08;
range = [2 253];
%% Optimization & Simulation
%for nr = 1:NR
    nr = 1;
    H = psrGrayHistogram(im);
    
    %Initial Krills positions
    for z1 = 1:NP
        X(z1,:) = sort(randperm(254, pop_size) + 1);
    end
    for n =1:size(X,2)
       X(:,n) = sort(X(:,n));
    end
    %plot(X(1,:),X(2,:),'o');
    for z2 = 1:pop_size
        K(z2)=cost(X(:,z2),H, q);
    end
    Kib=K;
    Xib=X;
    [Kgb(1,nr), A] = max(K);
    Xgb(:,1,nr) = X(:,A);
    p=1;  
    for j = 1:MI
        Xgb = sort(Xgb); 
        % Virtual Food
        for ll = 1:NP;
            if K ~= 0
                Sf(ll) = (sum(X(ll,:)./K));
            else
                Sf(ll) = 0;
            end
        end
        Xf(:,j) = Sf./(sum(1./K)); %Food Location     

        Xf(:,j) = findlimits(Xf(:,j)',LB,UB,Xgb(:,j,nr)'); %Bounds Checking
        Kf(j) = cost(Xf(:,j),H, q);
        if 2<=j
            if Kf(j-1)<Kf(j)
                Xf(:,j) = Xf(:,j-1);
                Kf(j) = Kf(j-1);
            end
        end
        
        Kw_Kgb = max(K)-Kgb(j,nr);
        % original w = (0.1+0.8*(1-j/MI));
        w = (0.1+0.8*(1-j/MI));
        
        for i = 1:pop_size
            % Calculation of distances
            Rf = Xf(:,j)-X(:,i);
            Rgb = Xgb(:,j,nr)-X(:,i);
            for ii = 1:pop_size
                RR(:,ii) = X(:,ii)-X(:,i);
            end
            R = sqrt(sum(RR.*RR));
                        
            % % % % % % % % % % % % % Movement Induced % % % % % % % % % %
            % Calculation of BEST KRILL effect
            if Kgb(j,nr) < K(i) && Kw_Kgb ~= 0
                alpha_b = -2*(1+rand*(j/MI))*(Kgb(j,nr) - K(i)) /Kw_Kgb/ sqrt(sum(Rgb.*Rgb)) * Rgb;
            else
                alpha_b=0;
            end
            
            %Calculation of NEIGHBORS KRILL effect
            nn=0;
            ds = mean(R)/5;
            alpha_n = 0;
            for n=1:pop_size
                if and(R<ds,n~=i)
                    nn=nn+1;
                    if and(nn<=4,K(i)~=K(n))
                        alpha_n = alpha_n-(K(n) - K(i)) /Kw_Kgb/ R(n) * RR(:,n);
                    end
                end
            end
            % Movement Induced
            N(:,i) = w*N(:,i)+Nmax*(alpha_b+alpha_n);
                        
            % % % % % % % % % % % % % Foraging Motion % % % % % % % % % %
            % Calculation of FOOD attraction
            if Kf(j) < K(i) && Kw_Kgb ~= 0 && sqrt(sum(Rf.*Rf)) ~= 0
                %'Kf(j)'
                %Kf(j)
                %'K(i)'
                %K(i)
                Beta_f = -2*(1-j/MI)*(Kf(j) - K(i))/Kw_Kgb/sqrt(sum(Rf.*Rf)) * Rf;
            else
                Beta_f=0;
            end
           
            %Beta_f = 0.5 * Beta_f; 
     
            
            % Calculation of BEST position attraction
            Rib = Xib(:,i)-X(:,i);
            if Kib(i) < K(i) && Kw_Kgb ~=0
                Beta_b=-(Kib(i) - K(i)) /Kw_Kgb/ sqrt(sum(Rib.*Rib)) *Rib;
            else
                Beta_b=0;
            end
            % Foraging Motion
            F(:,i) = w*F(:,i)+Vf*(Beta_b+Beta_f);
                        
            % % % % % % % % % % % % % Physical Diffusion % % % % % % % % %
            D = Dmax*(1-j/MI)*floor(rand+(K(i)-Kgb(j,nr))/Kw_Kgb)*(2*rand(NP,1)-ones(NP,1));
                        
            % % % % % % % % % % % % % Motion Process % % % % % % % % % % %

            DX = Dt*(N(:,i) + F(:,i));  
            %i
            %-(Kib(i) - K(i))
            %sqrt(sum(Rib.*Rib))
            %Rib
            %sqrt(sum(Rf.*Rf))
            %Rf
            %sqrt(sum(Rgb.*Rgb)) 
            %Rgb
            %DX
            %Beta_b
            %Beta_f
            %alpha_b
            %alpha_n
            % % % % % % % % % % % % % Crossover % % % % % % % % % % % % %
            if C_flag ==1
                C_rate = 0.8 + 0.2*(K(i)-Kgb(j,nr))/Kw_Kgb;
                Cr = rand(NP,1) < C_rate ;
                % Random selection of Krill No. for Crossover
                pop_size4Cr = round(pop_size*rand+.5);  
                % Crossover scheme
                X(:,i)=X(:,pop_size4Cr).*(1-Cr)+X(:,i).*Cr;
            end
            % Update the position
            X(:,i)=X(:,i)+DX;
            X(:,i)=sort(findlimits(X(:,i)',LB,UB,Xgb(:,j,nr)')); % Bounds Checking
              % X(:,i)  
            K(i)=cost(X(:,i),H, q);
            if K(i) > Kib(i)
                Kib(i)=K(i);
                Xib(:,i)=X(:,i);
            end     
                
        end
                
        % Update the current best
        [Kgb(j+1,nr), A] = max(K);
        if Kgb(j+1,nr)>Kgb(j,nr)
            Xgb(:,j+1,nr) = X(:,A);
        else
            Kgb(j+1,nr) = Kgb(j,nr);
            Xgb(:,j+1,nr) = Xgb(:,j,nr);
        end
        
        
        generation_entropy(j, 1) = Kgb(j,1);
        if j > 30  && std(generation_entropy(j-30:j))  < 0.01
           num_generations = j;
           break;
        end

    end
  
    %figure;
    %plot(X(1,:),X(2,:),'o');
    %hold on;
    %Sf = [0.05 0.05];
    %plot(Sf(:,1),Sf(:,2),'r*');
    %axis([xmin xmax ymin ymax]);
    


%% Post-Processing
[Best, Ron_No] = max(Kgb(end,:));
Lim = round(Xgb(:,end,Ron_No));
entropy = Kgb(end, 1);
%Mean = mean(Kgb(end,:))
%Worst = max(Kgb(end,:))
%Standard_Deviation = std(Kgb(end,:))

% Convergence plot of the best run
%figure;
%semilogy(1:MI+1,Kgb(:,Ron_No),1:MI+1,mean(Kgb'))
%xlabel('{\itNo. of Iterations}')
%ylabel('{\itf}({\bfx_{best}})')
%legend('Best run values','Average run values')

%figure;
%plot(X(1,:),X(2,:),'o');



function [ns]=findlimits(ns,Lb,Ub,best)
% Evolutionary Boundary Constraint Handling Scheme
n=size(ns,1);
for i=1:n
    ns_tmp=ns(i,:);
    I=ns_tmp<Lb;
    J=ns_tmp>Ub;
    A=rand;
    ns_tmp(I)=A*Lb(I)+(1-A)*best(I);
    B=rand;
    ns_tmp(J)=B*Ub(J)+(1-B)*best(J);
    ns(i,:)=ns_tmp;
end