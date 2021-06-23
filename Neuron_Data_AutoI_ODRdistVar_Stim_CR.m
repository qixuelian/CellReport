clear
[Neurons_num Neurons_txt] = xlsread('C:\work\DataBase\StimulationFilename_neuron_ODRdistVar.xlsx','neuronNoSigDecreaseEnoughTrial');%'BlockWise5_all''temp2'sigdecreaseAnyNoMixEnough  allneurons
Neurons = [Neurons_txt(:,1) num2cell(Neurons_num(:,1))];
IsMUA = 0;
phase='cuerate';
for n = 1:length(Neurons)
    if IsMUA
        filename = [Neurons{n,1},'_ch',num2str(Neurons{n,2})];
    else
        filename = [Neurons{n,1},'_',num2str(Neurons{n,2})];% filenameErr =[filename '_err'];
    end
%     Neuron_Data_Rasters_ODRdistVar_Cue_Stim_lineplot_2(filename);%plot
    
% %     phase =Neurons_num(n,2:3) ;
    [n_IndexC(n,:) fixrate(n,:)] = Neuron_Data_INDEX_odrdistVar_Stim(filename);
%     [n_IndexS(n,:) fixrate(n,:)] = Neuron_Data_INDEX_odrdistVar_sample_Stim(filename);
%     [n_IndexC200(n,:) fixrate(n,:)] = Neuron_Data_INDEX_odrdistVar_Stim200(filename);
%      [n_Index(n,:)  fixrate(n,:)] = Neuron_Data_INDEXofTask2_odrdistVar_Stim(filename);
%     [n_Index(n,:)] = Neuron_Data_INDEXofTask_odrdistVar_Stim200(filename,phase);
end
disp('done')
