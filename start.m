%
% This algorithm was written by Ricardo Morello Santos, Guilherme Wachs Lopes and
% Paulo Sérgio Rodrigues for his research in the area of image segmentation. 
% You can freely use it, modify it or pass forward for those you wish, since 
% it is intendend for using for academic and research purposes. 
%
% Bruteforce algorithm for testing bio-inspired algorithms in the BRATS
% database. We used the 2017 and 2018 versions of the BRATS database.
%
% Authors: Ricardo Morello Santos, Prof. Guilherme Wachs Lopes, 
% Prof. Nilson Saito and Prof. Paulo Sérgio Rodrigues
% Institution: Group of Signal Processing, Centro Universitario da FEI,
% Sao Bernardo do Campo, Brazil
% contact: unifrsantos@fei.edu.br
% Date: 2019
% Things you should consider changing for your computer's path:
% Path of the optimized slices file (line 49)
% Dataset path (line 40)
%This is the main function
function start
    %For parallelizating the process
    if matlabpool('size')==0  
        matlabpool(4)
    end
    %Number of algorithms compared
    %In this case, 7 algorithms were compared: 
    %WOA - Whale Optimization Algorithm
    %FF - Firefly Algorithm
    %GOA - Grasshopper Optimization Algorithm
    %KH - Krill Herd Algorithm
    %EHO - Elephant Herding Optimization
    %GWO - Grey Wolf Optimizer
    %CS - Cuckoo Search via Lévy Flights
    num_algorithms = 7;
    %Number of generations 
    generations = 200;
    %Size of the sample used (number of individuals from BRATS considered)
    sample_size = 100;
    %Path for the BRATS files
    dataset_path = 'E:\BRATS_Original\imagens_com_GT\amostra';
    %Number of folders for the BRATS database
    %Each BRATS year is divided into High-grade and Low-Grade glioma
    %In this case, four sub_folders were considered
    sub_folders = dir2(dataset_path);
    %Algorithms used
    algorithms = {'CS','WOA','KH', 'EHO', 'GOA', 'GWO', 'FF'};
    %Path for the file containing the optimized slices for each individual
    %Change this for the path in your computer
    [id_patient, slice1, slice2] = textread('E:\ideais.txt','%s %d %d');
    %Progressbar for user feedback
    progressbar('Progresso Geral')
    %Global variables that will be used for the progressbar
    global t1
    global t2
    global it
    it = 1;
    %For each folder (2017-HGG, 2017-LGG, 2018-HGG and 2018-LGG)
    for count_folders=1:length(sub_folders)
        %Retrieving these folders' path
        patients_path = strcat(dataset_path, strcat('\', sub_folders(count_folders).name));
        %Retrieving the patients paths
        patients = dir2(patients_path);
        %For each patient inside the folder
        for count_patients=1:length(patients)
           t1 = length(patients);
           t2 = 1;
           %Progress percentage
           percentage = count_patients/(length(patients));
           optimized_slice1 = 0;
           optimized_slice2 = 0;
           %Getting the optimized slices for the patient
           for i=1:sample_size
               if  strcmp(id_patient(i), char(patients(count_patients).name))
                  optimized_slice1 = slice1(i);
                  optimized_slice2 = slice2(i);
               end
           end 
           %Retrieving the patients names
           individuo = strcat(patients_path, strcat( '\', patients(count_patients).name));
           %For each algorithm 
           for i=1:num_algorithms      
               %Retrieving the algorithm to be evaluated
               algorithm_function = algorithms(i);
               %User feedback
               [algorithm_function percentage]
               %Getting the algorithms parameters
               parameters = retorna_parametros(algorithm_function);
               %Evaluating the algorithm function
               feval('execute_algorithm', algorithm_function, generations, parameters, individuo, optimized_slice1, optimized_slice2);
           end 
        end 
    end
end