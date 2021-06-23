%from Neuron_Data_SlidingAnova_ODRdistVr_stim_nostim_task
clear
excelName = 'C:\work\DataBase\StimulationFilename_neuron_ODRdistVar.xlsx';
sheetName1 = 'neuronNoSigDecreaseEnoughTrial' %'allNeurons'%'sigIncreaseAnyNoMix'% 'allNeurons';
[Neurons_num Neurons_txt] = xlsread(excelName,sheetName1);
Neurons = [Neurons_txt(:,1) num2cell(Neurons_num(:,1))];

TI = 0.01; % sliding time interval(sec) for ROC
TW = 0.500 % Time Window (sec) for ROC 0.25
startt=-1.0
endt= 5
cuetime = startt:TI:endt;

stim =1;
nUse = 0;
for n = 1: length(Neurons)
    filename = [Neurons{n,1},'_',num2str(Neurons{n,2})];
    load(filename);
    
    if ~isempty(MatData)
        Nclass = length(MatData.class);
        nn=0;
        
        %check if enough stim trials in matdata
        index1 =[];classuse =[1 4 6 9];
        for nt= 1:length(classuse)
            index1(nt) = length(find([MatData.class(classuse(nt)).ntr.Stim]==1));
            index3(nt) = length(find([MatData.class(classuse(nt)).ntr.Stim]==0));
        end
        
        index2 =[];classuse =[11 14 16 19];
        for nt= 1:length(classuse)
            index2(nt) = length(find([MatData.class(classuse(nt)).ntr.Stim]==1));
            index4(nt) = length(find([MatData.class(classuse(nt)).ntr.Stim]==0));
        end
        
        % collectrate for 2-way anova, stimuli only shown in RF or opposite of RF
        if isempty(find(index1<2)) & isempty(find(index2<2)) & isempty(find(index3<2)) & isempty(find(index4<2))
            
            for t =(startt-TW/2)/TI:1:(endt-TW/2)/TI % t = (-1-TW/2)/TI:1:(1.5-TW/2)/TI
                gs1=[];gs2=[];gs3=[];
                ntrW=[];
                rateAll = [];
                gs1_stim=[];gs2_stim=[];gs3_stim=[];
                ntrW_stim=[];
                rateAll_stim = [];
                nn=nn+1;
                classuse =[1 4 6 9 11 14 16 19];
                for nc = 1:length(classuse)
                    classN = classuse(nc);
                    index =find([MatData.class(classN).ntr(:).Stim] == 0);
                    MatData1.class(classN).ntr = MatData.class(classN).ntr(index);
                    index =find([MatData.class(classN).ntr(:).Stim] == 1);
                    MatData1_stim.class(classN).ntr = MatData.class(classN).ntr(index);%no stim trials
                    for i = 1:length(MatData1.class(classN).ntr)
                        Cue_onT1 = MatData1.class(classN).ntr(i).Cue_onT;
                        ntr_TS1 = MatData1.class(classN).ntr(i).TS ;
                        TS1 = ntr_TS1(find(ntr_TS1 >= Cue_onT1+(TI*t) & ntr_TS1 < Cue_onT1+TI*t+TW));
                        rateAll =  [rateAll length(TS1)/TW];
                        if classN <5 | (classN >10 & classN<15)
                            gs1=[gs1 1]; % location of 1st sti: in RF
                        else
                            gs1=[gs1 0]; % location of 1st sti:  out RF
                        end
                        
                        if classN == 4 | classN == 6  | classN ==14 | classN == 16
                            gs2=[gs2 1]; % location of 2nd stimuli in RF
                        else
                            gs2=[gs2 0]; % location of 2nd stimuli out RF
                        end
                        if classN <11
                            gs3=[gs3 1]; % white
                        else
                            gs3=[gs3 0]; % blue
                        end
                        
                    end
                    ntrW=[ntrW i];
                    
                    %
                    for i = 1:length(MatData1_stim.class(classN).ntr)
                        Cue_onT1 = MatData1_stim.class(classN).ntr(i).Cue_onT;
                        ntr_TS1 = MatData1_stim.class(classN).ntr(i).TS ;
                        TS1 = ntr_TS1(find(ntr_TS1 >= Cue_onT1+(TI*t) & ntr_TS1 < Cue_onT1+TI*t+TW));
                        rateAll_stim =  [rateAll_stim length(TS1)/TW];
                        if classN <5 | (classN >10 & classN<15)
                            gs1_stim=[gs1_stim 1]; % location of 1st sti: in RF
                        else
                            gs1_stim=[gs1_stim 0]; % location of 1st sti:  out RF
                        end
                        
                        if classN == 4 | classN == 6  | classN ==14 | classN == 16
                            gs2_stim=[gs2_stim 1]; % location of 2nd stimuli in RF
                        else
                            gs2_stim=[gs2_stim 0]; % location of 2nd stimuli out RF
                        end
                        
                        if classN <11
                            gs3_stim=[gs3_stim 1]; % white
                        else
                            gs3_stim=[gs3_stim 0]; % blue
                        end
                    end
                    ntrW_stim=[ntrW_stim i];
                end
                
                [pT tb stats]= anovan( rateAll,[gs1' gs2' gs3'], 'model', 'full','display','off');
                %              [pT tb stats]= anova1( rateAll,[gs1'], 'display','off');
                pW1(n,nn)=pT(1);
                pW2(n,nn)=pT(2);
                pW3(n,nn)=pT(3);
                pW4(n,nn)=pT(4);
                pW5(n,nn)=pT(5);
                pW6(n,nn)=pT(6);
                PEV1(n,nn) =  tb{2,2}/tb{10,2};%*100./length(rateAll);
                PEV2(n,nn) =  tb{3,2}/tb{10,2};%*100./length(rateAll);
                PEV3(n,nn) =  tb{4,2}/tb{10,2};%*100./length(rateAll);
                %                 SumSqW1(n,nn) = tb{2,2}/tb{6,2}*100./length(rateAll);
                %                 SumSqW2(n,nn) = tb{3,2}/tb{6,2}*100./length(rateAll);
                %                 SumSqW3(n,nn) = tb{4,2}/tb{6,2}*100./length(rateAll);
                
                [pT tb stats]= anovan( rateAll_stim,[gs1_stim' gs2_stim' gs3_stim' ], 'model', 'full','display','off');
                %              [pT tb stats]= anova1( rateAll,[gs1'], 'display','off');
                pW1_stim(n,nn)=pT(1);
                pW2_stim(n,nn)=pT(2);
                pW3_stim(n,nn)=pT(3);
                pW4_stim(n,nn)=pT(4);
                pW5_stim(n,nn)=pT(5);
                pW6_stim(n,nn)=pT(6);
                PEV1_stim(n,nn) =  tb{2,2}/tb{10,2};%*100./length(rateAll_stim);
                PEV2_stim(n,nn) =  tb{3,2}/tb{10,2};%*100./length(rateAll_stim);
                PEV3_stim(n,nn) =  tb{4,2}/tb{10,2};%*100./length(rateAll_stim);
                %                 SumSqW1_stim(n,nn) = tb{2,2}/tb{6,2}*100./length(rateAll_stim);
                %                 SumSqW2_stim(n,nn) = tb{3,2}/tb{6,2}*100./length(rateAll_stim);
                %                 SumSqW3_stim(n,nn) = tb{4,2}/tb{6,2}*100./length(rateAll_stim);
            end
            AUCss1_stim(n,1:length(cuetime))=smooth(pW1_stim(n,:),5);
            AUCss2_stim(n,1:length(cuetime))=smooth(pW2_stim(n,:),5);
            AUCss3_stim(n,1:length(cuetime))=smooth(pW3_stim(n,:),5)';
            AUCss1(n,1:length(cuetime))=smooth(pW1(n,:),5)';
            AUCss2(n,1:length(cuetime))=smooth(pW2(n,:),5)';
            AUCss3(n,1:length(cuetime))=smooth(pW3(n,:),5)';
            nUse = nUse+1;
            filenameUse{nUse,1}=Neurons{n,1};
            filenameUse{nUse,2}=     Neurons{n,2};
        else
            pW1(n,1:length(cuetime))=nan;
            pW2(n,1:length(cuetime))=nan;
            pW3(n,1:length(cuetime))=nan;
            pW4(n,1:length(cuetime))=nan;
            pW5(n,1:length(cuetime))=nan;
            pW6(n,1:length(cuetime))=nan;
            PEV1(n,1:length(cuetime)) =  nan;
            PEV2(n,1:length(cuetime)) =  nan;
            PEV3(n,1:length(cuetime)) =nan;
            
            %             SumSqW1(n,1:125) =nan;
            %             SumSqW2(n,1:125) = nan;
            %             SumSqW3(n,1:125) = nan;
            %             pB1(n,1:125)=nan;
            %             pB2(n,1:125)=nan;
            %             pB3(n,1:125)=nan;
            %             SumSqB1(n,1:125) =nan;
            %             SumSqB2(n,1:125) = nan;
            %             SumSqB3(n,1:125) = nan;
            %
            pW1_stim(n,1:length(cuetime))=nan;
            pW2_stim(n,1:length(cuetime))=nan;
            pW3_stim(n,1:length(cuetime))=nan;
            pW4_stim(n,1:length(cuetime))=nan;
            pW5_stim(n,1:length(cuetime))=nan;
            pW6_stim(n,1:length(cuetime))=nan;
            PEV1_stim(n,1:length(cuetime)) =  nan;
            PEV2_stim(n,1:length(cuetime)) =  nan;
            PEV3_stim(n,1:length(cuetime)) =nan;
            
            AUCss1_stim(n,1:length(cuetime))=nan;
            AUCss2_stim(n,1:length(cuetime))=nan;
            AUCss3_stim(n,1:length(cuetime))=nan;
            AUCss1(n,1:length(cuetime))=nan;
            AUCss2(n,1:length(cuetime))=nan;
            AUCss3(n,1:length(cuetime))=nan;
            
            %             SumSqW1_stim(n,1:125) =nan;
            %             SumSqW2_stim(n,1:125) = nan;
            %             SumSqW3_stim(n,1:125) = nan;
            %             pB1_stim(n,1:125)=nan;
            %             pB2_stim(n,1:125)=nan;
            %             pB3_stim(n,1:125)=nan;
            %             SumSqB1_stim(n,1:125) =nan;
            %             SumSqB2_stim(n,1:125) = nan;
            %             SumSqB3_stim(n,1:125) = nan;
        end
    end
end
% [a b c ] = find(pW1<0.05);

for n = 1 : length(cuetime)
    ntrpW1=  length(find(~isnan(pW1(:,1))));
    %     ntrpB1=  length(find(~isnan(pB1(:,1))));
    NpW1(n)=length(find(pW1(:,n)<0.05))./ntrpW1;
    NpW2(n)=length(find(pW2(:,n)<0.05))./ntrpW1;
    NpW3(n)=length(find(pW3(:,n)<0.05))./ntrpW1;
    NpW4(n)=length(find(pW4(:,n)<0.05))./ntrpW1;
    NpW5(n)=length(find(pW5(:,n)<0.05))./ntrpW1;
    NpW6(n)=length(find(pW6(:,n)<0.05))./ntrpW1;
    %     NpB1(n)=length(find(pB1(:,n)<0.05))./ntrpB1;
    %     NpB2(n)=length(find(pB2(:,n)<0.05))./ntrpB1;
    %     NpB3(n)=length(find(pB3(:,n)<0.05))./ntrpB1;
    
    
    ntrpW1_stim=  length(find(~isnan(pW1_stim(:,1))));
    %     ntrpB1_stim=  length(find(~isnan(pB1_stim(:,1))));
    NpW1_stim(n)=length(find(pW1_stim(:,n)<0.05))./ntrpW1_stim;
    NpW2_stim(n)=length(find(pW2_stim(:,n)<0.05))./ntrpW1_stim;
    NpW3_stim(n)=length(find(pW3_stim(:,n)<0.05))./ntrpW1_stim;
    NpW4_stim(n)=length(find(pW4_stim(:,n)<0.05))./ntrpW1_stim;
    NpW5_stim(n)=length(find(pW5_stim(:,n)<0.05))./ntrpW1_stim;
    NpW6_stim(n)=length(find(pW6_stim(:,n)<0.05))./ntrpW1_stim;
    %     NpB1_stim(n)=length(find(pB1_stim(:,n)<0.05))./ntrpB1_stim;
    %     NpB2_stim(n)=length(find(pB2_stim(:,n)<0.05))./ntrpB1_stim;
    %     NpB3_stim(n)=length(find(pB3_stim(:,n)<0.05))./ntrpB1_stim;
    
end
% pSumSqw1=nanmean(SumSqW1,1);
% pSumSqw2=nanmean(SumSqW2,1);
% pSumSqw3=nanmean(SumSqW3,1);
% pSumSqb1=nanmean(SumSqB1,1);
% pSumSqb2=nanmean(SumSqB2,1);
% pSumSqb3=nanmean(SumSqB3,1);

pSumSq1=nanmean(PEV1,1);
pSumSq2=nanmean(PEV2,1);
pSumSq3=nanmean(PEV3,1);
pSumSq1_stim=nanmean(PEV1_stim,1);
pSumSq2_stim=nanmean(PEV2_stim,1);
pSumSq3_stim=nanmean(PEV3_stim,1);
% pSumSq1=pSumSq1(1:121);
% pSumSq2=pSumSq2(1:121);
% pSumSq3=pSumSq3(1:121);
% pSumSq1_stim=pSumSq1_stim(1:121);
% pSumSq2_stim=pSumSq2_stim(1:121);
% pSumSq3_stim=pSumSq3_stim(1:121);
%Stim
cue = 0;
sample = 1.5;
target = 3;

% figure;hold on;
% plot(cuetime,NpW1,'-k','LineWidth',2)
% plot(cuetime,NpW1_stim,'-r','LineWidth',2)
% legend({'nostim_white' 'stim_white'});
% line([cue cue],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
% line([cue+.5 cue+.5],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
% line([sample sample],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
% line([sample+.5 sample+.5],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
% line([target target],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
% title([sheetName1 ' Number of neuron 3-wayAnova sig loc 1']);
% axis([-1 5 0 0.7]);
% 
% 
% figure;hold on;
% plot(cuetime,NpW2,'-b','LineWidth',2)
% plot(cuetime,NpW2_stim,'-r','LineWidth',2)
% legend({'nostim_blue' 'stim_blue'});
% line([cue cue],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
% line([cue+.5 cue+.5],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
% line([sample sample],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
% line([sample+.5 sample+.5],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
% line([target target],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
% title([sheetName1 ' Number of neuron 3-wayAnova sig loc 2']);
% axis([-1 5 0 0.7]);
% 
% 
% figure;hold on;
% plot(cuetime,NpW3,'-k','LineWidth',2)
% plot(cuetime,NpW3_stim,'-r','LineWidth',2)
% legend({'nostim_white' 'stim_white'});
% line([cue cue],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
% line([cue+.5 cue+.5],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
% line([sample sample],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
% line([sample+.5 sample+.5],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
% line([target target],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
% title([sheetName1 ' Number of neuron 3-wayAnova sig task']);
% axis([-1 5 0 0.7]);


figure;hold on;
plot(cuetime,pSumSq1,'-k','LineWidth',2)
plot(cuetime,pSumSq1_stim,'-r','LineWidth',2)
legend({'nostim_white' 'stim_white'});
line([cue cue],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
line([cue+.5 cue+.5],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
line([sample sample],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
line([sample+.5 sample+.5],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
line([target target],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
title([sheetName1 ' PEV 3-wayAnova sig loc 1']);
axis([-1 5 0 0.7]);


figure;hold on;
plot(cuetime,pSumSq2,'-b','LineWidth',2)
plot(cuetime,pSumSq2_stim,'-r','LineWidth',2)
legend({'nostim_blue' 'stim_blue'});
line([cue cue],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
line([cue+.5 cue+.5],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
line([sample sample],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
line([sample+.5 sample+.5],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
line([target target],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
title([sheetName1 ' PEV 3-wayAnova sig loc 2']);
axis([-1 5 0 0.7]);


figure;hold on;
plot(cuetime,pSumSq3,'-k','LineWidth',2)
plot(cuetime,pSumSq3_stim,'-r','LineWidth',2)
legend({'nostim_white' 'stim_white'});
line([cue cue],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
line([cue+.5 cue+.5],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
line([sample sample],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
line([sample+.5 sample+.5],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
line([target target],[0 100],'linestyle','--','LineWidth', 2, 'Color', 'k')
title([sheetName1 ' PEV 3-wayAnova sig task']);
axis([-1 5 0 0.7]);



