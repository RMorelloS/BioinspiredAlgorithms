% SEGMENTATION USING CS MCCULLOCH  ALGO WITH KAPUR'S ENTROPY AS OBJ FUNC: FOR RGB/GRAY IMAGES
%
% -----------------------------------------------------------------
% This demo program implements a modified version of Cuckoo Search (CS),   %
% using McCulloch's algorithm for Levy flights  generation of new solutions%
% coded by Shilpa Suresh.
% This is a modification of  the standard Cuckoo Search (CS) algorithm     %
% by Xin-She Yang at Cambridge University.
% -----------------------------------------------------------------
% Papers -- Citation Details:
% 1) Suresh, Shilpa, and Shyam Lal. "An efficient cuckoo search algorithm
%    based multilevel thresholding for segmentation of satellite images 
%    using different objective functions."
%    Expert Systems with Applications  58 (2016): 184-209.

% 2) X.-S. Yang, S. Deb, Cuckoo search via Levy flights,
%    in: Proc. of World Congress on Nature & Biologically Inspired
%    Computing (NaBIC 2009), December 2009, India,
%    IEEE Publications, USA,  pp. 210-214 (2009).
%
% 3) X.-S. Yang, S. Deb, Engineering optimization by cuckoo search,
%    Int. J. Mathematical Modelling and Numerical Optimisation, 
%    Vol. 1, No. 4, 330-343 (2010). 
%
% --------------------------------------------------------------- %


function [bestnestR,fmax] = CS2(histograma, thresholds, geracoes, q, parametros)
tic;
%if nargin<1,
% Number of nests (or different solutions)
n=parametros.n;%(i.e cuckoos( new solution) can lay eggs in any of these n nest)
%end
 
% Discovery rate of alien eggs/solutions
pa=parametros.pa;%(how well the host birdscan detect alian eggs)
 
% Change this if you want to get better results
N_IterTotalR=geracoes;
N_IterTotalG=geracoes;
N_IterTotalB=geracoes; 
%Data
%I=imread('image.jpg');
% I=rgb2gray(I);
Lmax= parametros.UB;
LMin = parametros.LB;
%Nt=size(I,1)*size(I,2);


% Simple bounds of the search domain
nd=thresholds;% number of thresholds required 

probR = psrGrayHistogram(histograma); 


%Lower and upper bounds
LbR=LMin*ones(1,nd); 
UbR=Lmax*ones(1,nd);  %(here it is from -5 to 5)
 
fitnessR=zeros(n,1);
% Random initial solutions
for i=1:n,
nestR(i,:)=sort(randperm(256,3)-1);% size(nest)=25X5....ie nXnd
end
for si=1:length(nestR)%No. of rows=N..ie similar si=1:N
nestR(si,:)=sort(nestR(si,:)); % sorting the xR generated randomly as above each row in ascending order
end

% initialized the population with n=25 nests, and nd=5 birds in each
% nest.the Lb and Ub decides the value bounds that can be assigned to each
% bird.
% Get the current best(finding the fittest one in each nest)
[fmaxR,bestnestR,nestR,fitnessR]=get_best_nest(nestR,nestR,fitnessR,nd,probR, q);


N_iterR=0;
N_iterG=0;
N_iterB=0;
% Starting iterations
for iter=1:N_IterTotalR,
    
     %Generate new solutions (but keep the current best)    
     %new_nestR=get_cuckoos(nestR,bestnestR,LbR,UbR,nd)

     new_nestR=empty_nests(nestR,LbR,UbR,pa);
     
     %[nestR new_nestR]
     %sort(sum(sqrt((nestR - new_nestR).^2)')')
     
     [fmax1R,bestR,nestR,fitnessR]=get_best_nest(nestR,new_nestR,fitnessR,nd,probR, q);
     
     T = [new_nestR fitnessR];     
     T = -sortrows(-T,4);
     new_nestR = T(:,1:3);
     fitnessR = T(:,4);
     
    %Update the counter(after one go through all nests)
     N_iterR=N_iterR+n; 
    %Discovery and randomization
    %Evaluate fitness for this set of solutions
    %[fmax1R,bestR,nestR,fitnessR]=get_best_nest(nestR,new_nestR,fitnessR,nd,probR, q);
    %Update the counter again
     % N_iterR=N_iterR+n;
    % Find the best objective so far  
    if fmax1R>fmaxR,
        fmaxR=fmax1R;
        bestnestR= bestR;
    end

end %% End of iterations


% Displaying segmented output
bestR=sort(bestR);
%Iout=imageGRAY(I,bestR);
bestnest=bestnestR;     %return optimal intensity
fmax=fmaxR;

time=toc; 
 
end



function imgOut=imageRGB(img,Rvec,Gvec,Bvec)%img=original image;Rvec=xR;Gvec=xG,Bvec=xB
imgOutR=img(:,:,1);
imgOutG=img(:,:,2);
imgOutB=img(:,:,3);
 
Rvec=[0 Rvec 256];
for iii=1:size(Rvec,2)-1
    at=find(imgOutR(:,:)>=Rvec(iii) & imgOutR(:,:)<Rvec(iii+1));
    imgOutR(at)=Rvec(iii);
end
 
Gvec=[0 Gvec 256];
for iii=1:size(Gvec,2)-1
    at=find(imgOutG(:,:)>=Gvec(iii) & imgOutG(:,:)<Gvec(iii+1));
    imgOutG(at)=Gvec(iii);
end
 
Bvec=[0 Bvec 256];
for iii=1:size(Bvec,2)-1
    at=find(imgOutB(:,:)>=Bvec(iii) & imgOutB(:,:)<Bvec(iii+1));
    imgOutB(at)=Bvec(iii);
end
 
imgOut=img;
 
imgOut(:,:,1)=imgOutR;
imgOut(:,:,2)=imgOutG;
imgOut(:,:,3)=imgOutB;
 
end

 function imgOut=imageGRAY(img,Rvec)
% imgOut=img;
limites=[0 Rvec 255];
tamanho=size(img);
imgOut(:,:)=img*0;
% cores=[ 0   0   0;
%         255 0   0;
%         0   255 0;
%         0   0   255;
%         255 255 0;
%         0   255 255;
%         255 0   255;
%         255 255 255];
        
cores=colormap(lines)*255;
close all;
%tic
k=1;
    for i= 1:tamanho(1,1)
        for j=1:tamanho(1,2)
            while(k<size(limites,2))
                if(img(i,j)>=limites(1,k) && img(i,j)<=limites(1,k+1))
                    imgOut(i,j,1)=limites(1,k);
%                     imgOut(i,j,2)=cores(k,2);
%                     imgOut(i,j,3)=cores(k,3);
                end
                k=k+1;
            end
            k=1;
        end
    end
 
 end
 
%% --------------- All subfunctions are list below ------------------
%% Get cuckoos by random walk
function nest=get_cuckoos(nest,best,Lb,Ub,nd)

% Levy flights
n=size(nest,1);
% Levy exponent and coefficient
% For details, see equation (2.21), Page 16 (chapter 2) of the book
% X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).
beta=3/2;

sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);

for j=1:n,
   
   s=nest(j,:);
   % This is a simple way of implementing Levy flights
   % For standard random walks, use step=1;
   
   %% Levy flights by Mcculloch's algorithm
   x = stabrnd(.5,1, 1, 1,1, nd);
   
   % Now the actual random walks or flights
   s=s+x;
   
   % Apply simple bounds/limits
   nest(j,:)=simplebounds(s,Lb,Ub);% assuring the new updated solution is within the bounds and replacing the corresponding nest values
   
   for si=1:n   %No. of rows=N..ie similar si=1:N
      nest(si,:)=sort(nest(si,:)); % sorting the xR generated randomly as above each row in ascending order
   end
   
   nest(j,:)=fix(nest(j,:));

end

end



%% Find the current best nest
function [fmax,best,nest,fitness]=get_best_nest(nest,newnest,fitness,nd,probR, q)
% Evaluating all new solutions

for j=1:size(nest,1),
    fnew=fobj(newnest(j,:),nd,probR, q);
    if fnew>=fitness(j),
       fitness(j)=fnew;
       nest(j,:)=newnest(j,:);
    end
end

% Find the current best
[fmax,K]=max(fitness); %fmin=minimum fitness;K=corresponding row.ie nest
best=nest(K,:);
end
 
%% Replace some nests by constructing new solutions/nests
function new_nest=empty_nests(nest,Lb,Ub,pa)
% A fraction of worse nests are discovered with a probability pa
n=size(nest,1);
% Discovered or not -- a status vector
K=rand(size(nest))>pa;
 
% In the real world, if a cuckoo's egg is very similar to a host's eggs, then 
% this cuckoo's egg is less likely to be discovered, thus the fitness should 
% be related to the difference in solutions.  Therefore, it is a good idea 
% to do a random walk in a biased way with some random step sizes.  
%% New solution by biased/selective random walks
stepsize=rand*(nest(randperm(n),:)-nest(randperm(n),:));
new_nest=nest+stepsize.*K;
for j=1:size(new_nest,1)
    s=new_nest(j,:);
    new_nest(j,:)=simplebounds(s,Lb,Ub); 
end    
for si=1:n%No. of rows=N..ie similar si=1:N
nest(si,:)=sort(new_nest(si,:)); % sorting the xR generated randomly as above each row in ascending order
end
new_nest=fix(nest);

end


% Application of simple constraints
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound
  ns_tmp=s;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);
  
  % Apply the upper bounds 
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move 
  s=ns_tmp;
 
%% fitness function
end


   
function fnew=fobj(u,nd,probR, q)

    if q == 1
        fnew = psrAvaliacaoShannon(probR, u);
    else
        fnew = psrAvaliacaoTsallis(probR, q, u);
    end
end
    
    
function light = shannon(u, nd, probR)
light = 0;
n = size(u,2) + 2;
elemento=[0 u 255];

a = elemento(1);
b = elemento(2);

light = ShannonEntropy(probR,a,b);

for i=2:n-1
  a = elemento(i) + 1;
  b = elemento(i+1);
  ES = ShannonEntropy(probR,a,b);
  light = light + ES;
end
end




function S = ShannonEntropy(X,a,b)
  a = round(a);
  b = round(b);
  if b == 255
      b = 254;
  end
  H = zeros(1,b-a+1);
  H = X(a + 1 :b + 1);
  H = H/sum(H);
  L = size(H,2);
  S = 0;
  ret = 0;
  for i=1:L
    if H(i) ~= 0
      ret = ret + H(i) * log(H(i));
      if isnan(ret)
         ret = 0;
      end
    end
  end 
  
  S = -ret;
end 



 % Stable Random Number Generator (McCulloch 12/18/96)
%------------------------------------------------------------------------------
%------------------------------------------------------------------------------

 function [x] = stabrnd(alpha, beta, c, delta, m, n)

 % Returns m x n matrix of iid stable random numbers with
 % characteristic exponent alpha in [.1,2], skewness parameter
 % beta in [-1,1], scale c > 0, and location parameter delta.
 % Based on the method of J.M. Chambers, C.L. Mallows and B.W.
 % Stuck, "A Method for Simulating Stable Random Variables,"
 % JASA 71 (1976): 340-4.
 % Encoded in MATLAB by J. Huston McCulloch, Ohio State
 % University Econ. Dept. (mcculloch.2@osu.edu). This 12/18/96
 % version uses 2*m*n calls to RAND, and does not rely on
 % the STATISTICS toolbox. 
 % The CMS method is applied in such a way that x will have the
 % log characteristic function
 % log E exp(ixt) = i*delta*t + psi(c*t),
 % where
 % psi(t) = -abs(t)^alpha*(1-i*beta*sign(t)*tan(pi*alpha/2))
 % for alpha ~= 1,
 % = -abs(t)*(1+i*beta*(2/pi)*sign(t)*log(abs(t))),
 % for alpha = 1.
 % With this parameterization, the stable cdf S(x; alpha, beta,
 % c, delta) equals S((x-delta)/c; alpha, beta, 1, 0). See my
 % "On the parametrization of the afocal stable distributions,"
 % _Bull. London Math. Soc._ 28 (1996): 651-55, for details.
 % When alpha = 2, the distribution is Gaussian with mean delta
 % and variance 2*c^2, and beta has no effect.
 
 % When alpha > 1, the mean is delta for all beta. When alpha
 % <= 1, the mean is undefined.
 % When beta = 0, the distribution is symmetrical and delta is
 % the median for all alpha. When alpha = 1 and beta = 0, the
 % distribution is Cauchy (arctangent) with median delta.
 % When the submitted alpha is > 2 or < .1, or beta is outside
 % [-1,1], an error message is generated and x is returned as a
 % matrix of NaNs.
 % Alpha < .1 is not allowed here because of the non-negligible
 % probability of overflows.

 % If you're only interested in the symmetric cases, you may just
 % set beta = 0 and skip the following considerations:
 % When beta > 0 (< 0), the distribution is skewed to the right
 % (left).
 % When alpha < 1, delta, as defined above, is the unique fractile
 % that is invariant under averaging of iid contributions. I
 % call such a fractile a "focus of stability." This, like the
 % mean, is a natural location parameter.
 % When alpha = 1, either every fractile is a focus of stability,
 % as in the beta = 0 Cauchy case, or else there is no focus of
 % stability at all, as is the case for beta ~=0. In the latter
 % cases, which I call "afocal," delta is just an arbitrary
 % fractile that has a simple relation to the c.f.
 % When alpha > 1 and beta > 0, med(x) must lie very far below
 % the mean as alpha approaches 1 from above. Furthermore, as
 % alpha approaches 1 from below, med(x) must lie very far above
 % the focus of stability when beta > 0. If beta ~= 0, there
 % is therefore a discontinuity in the distribution as a function
 % of alpha as alpha passes 1, when delta is held constant.
 % CMS, following an insight of Vladimir Zolotarev, remove this
 % discontinuity by subtracting
 % beta*c*tan(pi*alpha/2)
 % (equivalent to their -tan(alpha*phi0)) from x for alpha ~=1
 % in their program RSTAB, a.k.a. RNSTA in IMSL (formerly GGSTA).
 % The result is a random number whose distribution is a contin-
 % uous function of alpha, but whose location parameter (which I
 % call zeta) is a shifted version of delta that has no known
 % interpretation other than computational convenience.
 % The present program restores the more meaningful "delta"
 % parameterization by using the CMS (4.1), but with
 % beta*c*tan(pi*alpha/2) added back in (ie with their initial
 % tan(alpha*phi0) deleted). RNSTA therefore gives different
 % results than the present program when beta ~= 0. However,
 % the present beta is equivalent to the CMS beta' (BPRIME).
 % Rather than using the CMS D2 and exp2 functions to compensate
 % for the ill-condition of the CMS (4.1) when alpha is very
 % near 1, the present program merely fudges these cases by
 % computing x from their (2.4) and adjusting for
 % beta*c*tan(pi*alpha/2) when alpha is within 1.e-8 of 1.
 % This should make no difference for simulation results with
 % samples of size less than approximately 10^8, and then
 % only when the desired alpha is within 1.e-8 of 1, but not
 % equal to 1.
 % The frequently used Gaussian and symmetric cases are coded
 % separately so as to speed up execution.

 % Additional references:
 % V.M. Zolotarev, _One Dimensional Stable Laws_, Amer. Math.
 % Soc., 1986.
 % G. Samorodnitsky and M.S. Taqqu, _Stable Non-Gaussian Random
 % Processes_, Chapman & Hill, 1994.
 % A. Janicki and A. Weron, _Simulaton and Chaotic Behavior of
 % Alpha-Stable Stochastic Processes_, Dekker, 1994.
 % J.H. McCulloch, "Financial Applications of Stable Distributons," 
 % _Handbook of Statistics_ Vol. 14, forthcoming early 1997.

 % Errortraps:
 if alpha < .1 | alpha > 2
 disp('Alpha must be in [.1,2] for function STABRND.')
 alpha
 x = NaN * zeros(m,n);
 return
 end
 if abs(beta) > 1
 disp('Beta must be in [-1,1] for function STABRND.')
 beta
 x = NaN * zeros(m,n);
 return
 end

% Generate exponential w and uniform phi:
 w = -log(rand(m,n));
phi = (rand(m,n)-.5)*pi;

 % Gaussian case (Box-Muller):
 if alpha == 2
 x = (2*sqrt(w) .* sin(phi));
 x = delta + c*x;
 return
 end

 % Symmetrical cases:
 if beta == 0
if alpha == 1 % Cauchy case
 x = tan(phi);
 else
 x = ((cos((1-alpha)*phi) ./ w) .^ (1/alpha - 1) ...
 .* sin(alpha * phi) ./ cos(phi) .^ (1/alpha));
 end

 % General cases:
 else
 cosphi = cos(phi);
 if abs(alpha-1) > 1.e-8
 zeta = beta * tan(pi*alpha/2);
 aphi = alpha * phi;
 a1phi = (1 - alpha) * phi;
 x = ((sin(aphi) + zeta * cos(aphi)) ./ cosphi) ...
 .* ((cos(a1phi) + zeta * sin(a1phi)) ...
 ./ (w .* cosphi)) .^ ((1-alpha)/alpha);
 else
 bphi = (pi/2) + beta * phi;
 x = (2/pi) * (bphi .* tan(phi) - beta * log((pi/2) * w ...
 .* cosphi ./ bphi));
 if alpha ~= 1
 x = x + beta * tan(pi * alpha/2);
 end
 end
end

 % Finale:
x = delta + c * x;
 return
% End of STABRND.M 
end