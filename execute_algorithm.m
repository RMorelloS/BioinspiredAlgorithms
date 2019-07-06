% Function for evaluating each bio-inspired algorithm 
% input parameters: Algorithm ID, number of generations, algorithms
% parameters, patient's path, first optimized slice, second optimized slice
% output: a '.txt' file containing the results from the bruteforce algorithm
% Things you should consider changing to your computer's path:
% Output path for the '.txt' file: Line 82


function execute_algorithm(funcao,geracoes, parametros, path, slice_otimizado1, slice_otimizado2)
%Global variables for the progressbar
global t1;
global t2;
global it;
%Reading the ground-truth file for the patient
[volumeGT, idIndividuo] = leituraGT(path);
%Reading the names of the four remaining files (T1w, T2w, T2-Flair and
%T1ce)
pastaArquivos = dir2(path);
%For each one of the four files mentioned above
for k=1:numel(pastaArquivos)
    %Name of the file
    nomeArquivo = pastaArquivos(k).name;
    %Array containing both optimized slices
    slicesEscolhidos = [slice_otimizado1,slice_otimizado2];
    %Reading the ground-truth file
    niiVolume = leituraNII(path, nomeArquivo);
    %For each optimized slice
    for i=1:size(slicesEscolhidos,2)
        %Ground-truth slice
        sliceGT = squeeze(volumeGT(:,:,slicesEscolhidos(i)))';
        sliceGT = double(sliceGT);
        %Normalizing the ground-truth to the [0,255] range
        if max(max(sliceGT)) - min(min(sliceGT)) ~= 0
            sliceGT = (sliceGT-min(min(sliceGT)))/(max(max(sliceGT)) - min(min(sliceGT)));
        else
            sliceGT = zeros(size(sliceGT));
        end
        sliceGT = round(sliceGT * 255);
        %Number of the thresholds to segment the slice
        thresholds = size(unique(sliceGT),1) - 1;
        %Reading the magnetic ressonance slice and normalizing it to the 
        %[0-255] range
        sliceNII = squeeze(niiVolume(:,:,slicesEscolhidos(i)))';
        sliceNII = double(sliceNII);

        if max(max(sliceNII)) - min(min(sliceNII)) ~=0
            sliceNII = (sliceNII-min(min(sliceNII)))/(max(max(sliceNII)) - min(min(sliceNII)));
        end
        sliceNII = round(255 * sliceNII);
        %For - Executing each algorithm 10 times for the algorithm's 
        %output variance measurement   
        for j=1:10
           %For - q values in the [0.1, 2] range
           parfor q=1:20
               %tic in order to store the algorithm's execution time
                tic
                if thresholds >= 1
                    %Executing the algorithm
                    [Lims, entropia, numGeracoes, melhores_avaliacoes] = feval(char(funcao),sliceNII, thresholds, geracoes, q/10, parametros);
                    %Rounding the thresholds to int
                    Lims = round(Lims);
                    %Calculating the metrics 
                    [dice, jaccard,dklr] = metricas(Lims, sliceNII, sliceGT);
                else
                    Lims = -1;
                    entropia = -1;
                    dice = -1;
                    jaccard = -1;
                    dklr = -1;
                end
                t = toc;
                %Storing the results in a struct
                resultado = preencher_struct(parametros,i,q,j,dice,jaccard,dklr,entropia,...
                                   nomeArquivo,t,thresholds, funcao, geracoes, Lims, idIndividuo, numGeracoes, melhores_avaliacoes);
                %Storing the struct with the results in an array
                array(q) = resultado;
           end
            %Updating the progressbar
            it = it+1;
            progressbar(it / (10*numel(pastaArquivos)*size(slicesEscolhidos, 2)*t1*t2 * 4 + 1));
            %Output file to store the results
            outputFile = fopen('C:\Users\Ricardo\Downloads\GOA_geral.txt', 'a+'); 
            %For - storing the struct array in the file
            for cont_array=1:size(array,2)
                fprintf(outputFile, '%d ', slicesEscolhidos(array(cont_array).slice));
                fprintf(outputFile, '%f ', array(cont_array).q/10);
                fprintf(outputFile, '%d ', array(cont_array).iteracaoAlg);
                fprintf(outputFile, '%f ', array(cont_array).dice);
                fprintf(outputFile, '%f ', array(cont_array).jaccard);
                fprintf(outputFile, '%f ', array(cont_array).dkl);
                fprintf(outputFile, '%f ', array(cont_array).entropia);
                fprintf(outputFile, '%s ', char(array(cont_array).nome_arquivo));
                fprintf(outputFile, '%f ', array(cont_array).tempo);
                fprintf(outputFile, '%d ', array(cont_array).num_limiares);
                fprintf(outputFile, '%s ', char(array(cont_array).algoritmo));
                fprintf(outputFile, '%d ', array(cont_array).num_geracoes);
                fprintf(outputFile, '%d ', array(cont_array).numGeracoes);
                if thresholds >= 1
                    for cont_lims=1:thresholds
                        fprintf(outputFile, '%d ', array(cont_array).limiares(cont_lims));   
                    end
                else
                     fprintf(outputFile, '-1 ');   
                end 
                fprintf(outputFile, '%s ', char(array(cont_array).id_paciente));
                value = struct2cell(array(cont_array).parametros);
                for cont_parametros=1:size(value,1)
                    valorParam = value{cont_parametros};
                    fprintf(outputFile,'%d ',valorParam);
                end
                
                fprintf(outputFile, '{ ');
                for iter=1:array(cont_array).numGeracoes
                    fprintf(outputFile, '%f ', array(cont_array).melhores_avaliacoes(iter, 1)); 
                end
                fprintf(outputFile, '}');
                fprintf(outputFile, '\r\n');
            end
            fclose(outputFile);
        end
    end
end

end
