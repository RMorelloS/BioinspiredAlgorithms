%
% This algorithm was written by Paulo Sergio Rodrigues for his research in 
% the area of ??image segmentation. You can freely use it, modify it or pass 
% forward for those you wish, since it is intendend for using for academic and 
% research purposes. We recommend that whenever use, please reference 
% the following papers or use the following latex format:
%
%@INPROCEEDINGS{rodrigues2009,
%  AUTHOR =       {P.S. Rodrigues and G.A. Giraldi},
%  TITLE =        {Computing the q-index for Tsallis nonextensive image segmentation},
%  BOOKTITLE =    {XXII Brazilian Symposium on Computer Graphics and Image Processing (SIBGRAPI)},
%  YEAR =         {2009},
%  pages =        {232-237},
%  address =      {Rio de Janeiro, RJ, Brazil},
%}
%
%@ARTICLE{rodrigues2011a,
%  AUTHOR =       {P. S. RODRIGUES and G. A. Giraldi},
%  TITLE =        {Improving the Non-Extensive Medical Image Segmentation Based on Tsallis Entropy},
%  JOURNAL =      {Pattern Analysis and Applications (Print)},
%  YEAR =         {2011},
%  volume =       {14},
%  pages =        {1-14},
%}
%
%@INBOOK{psergiobook2008,
%  AUTHOR =       {P. S. Rodrigues and G. A. Giraldi and J. Suri and  R. F. Chang},
%  editor =       {J. Suri and C. Kathuria and R. F. Chang and F. Molinari and A. Fenster},
%  TITLE =        {Automatic Classification of Breast Lesions in 3D Ultrasound Images},
%  series =       {{Advances in Diagnostic and Therapeutic Ultrasound Imaging}},
%  PUBLISHER =    {Artech House},
%  pages     =    {189-223},
%  chapter  =     {8},
%  year     =     {2008},
%  address  =     {Boston and London}
%}
%
% ALGORITHM FOR Computing GrayScale Histogram Image according to
% the following paper from IEEE data base:
% 
% P. S. Rodrigues and G. A. Giraldi. "Improving the non-extensive medical image
% segmentation based on Tsallis entropy". Pattern Analysis and Applications, 
% v. 14, p. 1-18, 2011.
%
% P. S. Rodrigues and G. A. Giraldi, "Computing the q-index for Tsallis 
% Nonextensive Image Segmentation. In: 2009 XXII Brazilian Symposium on 
% Computer Graphics and Image Processing (SIBGRAPI), 2009, Rio de Janiero.
% 2009 XXII Brazilian Symposium on Computer Graphics and Image Processing. p. 232.
%
% P. S. Rodrigues, G. A. Giraldi et all. "Automatic Classification of Breast 
% Lesions in 3D Ultrasound Images. In: J. Suri, C. Kathuria, R. F. Chang, 
% F. Molinari, A. Fenster. (Org.). In Book: Advances in Diagnostic and 
% Therapeutic Ultrasound Imaging. Boston and London: Artech House, 
% 2008, v. , p. 189-223.
% 
% INPUT PARAMETERS: img -> an grayscale image
%
% OUTPUT PARAMETER: a Histogram distribution H
%
% Author: Prof. Paulo Sergio Rodrigues
% Institution: Group of Signal Processing, Centro Universitario da FEI,
% Sao Paulo, Brazil
% contact: pslucano@gmail.com
% Date: 2010
function H = psrGrayHistogram(im)

[Lin, Col] = size(im);
H = zeros(1,256);
for y=1:Lin
    for x=1:Col
      H(im(y,x)+1) = H(im(y,x)+1) + 1;
    end
end
% normalization
H = H/(sum(H));





