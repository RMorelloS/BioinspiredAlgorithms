# BioInspiredAlgorithmsComparison
Comparison of bio-inspired algorithms for medical image segmentation using Tsallis entropy on the BRATS dataset (years 2017 and 2018) using MATLAB.
# Authors
Ricardo Morello Santos, Prof. Guilherme Wachs Lopes, Prof. Nilson Saito and Prof. Paulo Sérgio Rodrigues
# General instructions
The main file for this experiment is "start.m". Each algorithm is executed through the "execute_algorithm.m" file. 
Also, each algorithm may be executed independently for any kind of image. Follow the instructions in the header of each algorithm file to get to know better about the parameters considered. A description of each algorithm and their respective parameters is suplied below.

# Algorithms and their parameters:
All algorithms expect four parameters as input: 
The image to be segmented: *I*
Number of segmentation thresholds: *thresholds*
Number of generations: *generations*
Parameter struct: *parameters*

The names of these parameters may vary for each algorithm. The *parameters* struct contains dinamic parameters which vary for each algorithm. In this study, we considered the following algorithms and their respective parameters:

# Cuckoo Search via Lévy Flights (CS): 

 Population size (pop_size) = 40;
 
 Probability of host bird discovering the cucoo egg (pa) = 0.5;
 
 Upper Bound for image thresholding (UB) = 253;
 
 Lower Bound for image thresholding (LB) = 2.
# Whale Optimization Algorithm (WOA): 

 Population size (pop_size) = 30;
 
 Upper Bound for image thresholding (UB) = 253;
 
 Lower Bound for image thresholding (LB) = 2.
 
 # Krill Herd Algorithm (KH):
 
 Population size (pop_size) = 40;
 
 Upper Bound for image thresholding (UB) = 253;
 
 Lower Bound for image thresholding (LB) = 2.
 
 # Elephant Herding Optimization (EHO):
 
 Population size (pop_size) = 100;
 
 Number of clans (numClan) = 5;
 
 Elitism (Keep) = 2;
 
 Alpha (alpha) = 0.5;
 
 Beta (beta) = 0.1;
 
 Upper Bound for image thresholding (UB) = 253;
 
 Lower Bound for image thresholding (LB) = 2.
 
 # Grasshopper Optimization Algorithm (GOA):
 
 Population size (pop_size) = 30;
 
 Upper Bound for image thresholding (UB) = 253;
 
 Lower Bound for image thresholding (LB) = 2.
 
 # Grey Wolf Optimizer (GWO):
 
 Population size (pop_size) = 30;
 
 Upper Bound for image thresholding (UB) = 253;
 
 Lower Bound for image thresholding (LB) = 2.
 
 # Firefly Algorithm (FA):
 
 Population size (pop_size) = 50;
 
 Upper Bound for image thresholding (UB) = 253;
 
 Lower Bound for image thresholding (LB) = 2.

