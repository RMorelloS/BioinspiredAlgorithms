%_________________________________________________________________________%
%  Grasshopper Optimization Algorithm (GOA) source codes demo V1.0        %
%                                                                         %
%  Developed in MATLAB R2016a                                             %
%                                                                         %
%  Author and programmer: Seyedali Mirjalili                              %
%                                                                         %
%         e-Mail: ali.mirjalili@gmail.com                                 %
%                 seyedali.mirjalili@griffithuni.edu.au                   %
%                                                                         %
%       Homepage: http://www.alimirjalili.com                             %
%                                                                         %
%  Main paper: S. Saremi, S. Mirjalili, A. Lewis                          %
%              Grasshopper Optimisation Algorithm: Theory and Application %
%               Advances in Engineering Software , in press,              %
%               DOI: http://dx.doi.org/10.1016/j.advengsoft.2017.01.004   %
%                                                                         %
%_________________________________________________________________________%
%
%This function randomly initializes the position of agents in the search space.
function [X]=initializationGOA(N,dim,up,down)
  for i=1:N
       %X(:,i)= randperm(254, N) + 1;
       X(i,:) = randperm(253, dim) + 2;
  end

end
