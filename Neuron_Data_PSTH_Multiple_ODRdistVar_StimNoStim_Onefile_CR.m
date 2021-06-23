%  for spatial only
%  from Neuron_Data_PSTH_PFCPPC.m
% 'Neuron_Data_Psth_Pre_post.m' can't plot lengend.
% updated on feb 28 2010, XQ
% filerange = 'Spatial'
% MaxOrOpp is using the max friting rate or the opposite/minum firing rate
% max vs worst location
clear
normalize = 0; % 1: normalize data 0, no normalization
smooth = 0; % 1: smooth data with gaussian, 0: no smooth
IsSmooth = 0;
Swindow = 2; % if smooth, smooth window
isSeperateMNM = 0; % seperate Match and NonMatch classes: 1; if not 0
MaxOROpp = 'max'; % collect data from max firing rate: 'max', or opposite: 'opp'

Npsth = 1;
excelName = 'C:\work\DataBase\StimulationFilename_neuron_ODRdistVar.xlsx'%ODRdistVar_filename.xlsx';
sheetName1 = 'neuronNoSigDecreaseEnoughTrial' % 'allNeurons_use' 
%'noSigNeurons_use' %'neuronNoSigDecreaseEnoughConly' %'neuronNoSigDecreaseEnoughCD' % 'neuronNoSigDecreaseEnoughCD' %'Neuron_Cueonly' %'neuronUseSigIncreaseFix' 
%'neuronNoSigDecreaseEnoughTrial' %'NoSigAnyChange'%'sigdecreaseAnyNoMixEnough'%'neuronsigDecreseenoughtrial'
%'Sheet2''neuronUseIncreaseFix1Trial'% 'noSigNeurons'%'neuronUseNoSigChangeFix'%'neuronUseSigIncreaseFix_use'
%'neuronUsedecreaseFix'%'neuronUseIncreaseFix'%'SigChangeEnoughTrials'%BlockWise5_HEC'%'BlockWise5_all'
%'SigALLPFCnostim';%'SigTaskPPCpreB'%'sigcuedelayPPC'
%'neuronUseIncreaseFix2Trial'-at least two trials each class
%neuronUseSigIncreaseFix
% sheetName2 = 'SigAll';
Area = '';
monkey = { 'hec'};
whichPhase = 'cuedelay' % using Cuerate/CueDelay/samplerate/sampledelay activity to selecte classes, max class
lengends = {'1stM','2ndM','1stOPP','2ndOPP'};para.legendsMNMOPP = {'1stM','2ndM','1stOPP','2ndOPP'};
filerange = ['odrdist']; %'Spatial' 'Spatial'
isFeature = [0 0];
colors = {'b','r','g','m'};
colors1 = {[0.678 0.922 1], [1 0.6 0.784]};
figname = [sheetName1 '_' whichPhase ' ' filerange];
classDef ={'max-180';'max-90';'max-45';'max-0';'max-null';'opp-180';'opp-90';'opp-45';'opp-0';'opp-null';...
    'max-180B';'max-90B';'max-45B';'max-0B';'max-nullB';'opp-180B';'opp-90B';'opp-45B';'opp-0B';'opp-nullB';};


outPutDir = ['C:\work\Data_Analysis\APM_Data\DA_Data\'];
warning off MATLAB:divideByZero
cd C:\work\Data_Analysis\APM_Data\DA_Data

oppClasslist = [5 6 7 8 1 2 3 4 5];
cue = 1;
sample = 2.5;
target = 4;
bin_width = 0.05;% normal 0.05;
bin_step=0.05;
bin_edges= 0:bin_width:5;
max_resp = 35;
fixrate = [];
para.colors = colors;
para.legends = lengends;
% para.axisNumber = [0 5.5 0 35];
if normalize
    para.axisNumber = [0 5 0 5];
else
    para.axisNumber = [0 5 0 60];
end
para.bins = bin_edges + 0.5 * bin_width;
para.times = [1 1.5;3 3.5 ;5 nan];
MaxClassAll=[];
nFile=0;
figure;hold on;
for Nsheet =  1 : Npsth
    eval(['[Neurons_num Neurons_txt] = xlsread(excelName,sheetName' num2str(Nsheet) ');'])
    Neurons = [Neurons_txt(:,1) num2cell(Neurons_num(:,1))];
    
    for n = 1: length(Neurons)
        filename = [Neurons{n,1},'_',num2str(Neurons{n,2})];
        if isFeature(Nsheet)
            task =Neurons{n,3};
        end
        % NM stim; M no Stim
        [psthM psthNM fixrateMean psthCT1 psthCT2 isReverse maxpsth(n,1) maxCpsth(n,1) ntrStimT ] = extractPSTH_OdrdistVar_StimOneFile(filename,bin_edges,bin_width,IsSmooth,whichPhase);
        %         [psthM psthNM fixrateMean psthCT1 psthCT2 isReverse maxpsth(n,1) maxCpsth(n,1) ntrStimT ] = extractPSTH_OdrdistVar_StimOneFile_enoughTrial(filename,bin_edges,bin_width,IsSmooth,whichPhase);
        %        [psthM psthNM fixrateMean psthCT1 psthCT2 isReverse maxpsth(n,1) maxCpsth(n,1)] =  extractPSTH_OdrdistVar_StimOneFile_2ndStimuli(filename,bin_edges,bin_width,IsSmooth,whichPhase);
        
        ntr_stim(n,:)=ntrStimT;
        if isempty(find(ntrStimT<3))
            %use
            nFile = nFile +1;
            neuron_use(nFile)=n;
            
        end
        if ~isempty(psthCT1(1).psth) & ~isempty(psthCT2(1).psth)
            
            for m = 1 : size(psthCT1,2)% {[1:4]; [6:9]; [11:14]; [16:19]}
                if normalize
                    eval(['psthCM' num2str(m) '(n,:)=psthCT1(m).psth./maxCpsth(n,1);'])%nonStim
                    if isempty(psthCT2(m).psth)
                        eval(['psthCNM' num2str(m) '(n,:)=nan;'])%stim
                    else
                        eval(['psthCNM' num2str(m) '(n,:)=psthCT2(m).psth./maxCpsth(n,1);'])%stim
                    end
                else
                    eval(['psthCM' num2str(m) '(n,:)=psthCT1(m).psth;'])%nonStim
                    if isempty(psthCT2(m).psth)
                        eval(['psthCNM' num2str(m) '(n,:)=nan;'])%stim
                    else
                        eval(['psthCNM' num2str(m) '(n,:)=psthCT2(m).psth;'])%stim
                    end
                    
                    % 2019 Apr 1
                    %                         max_resp = max([psthCT1(m).psth psthCT2(m).psth]);
                    %                         if ~isnan(max_resp)
                    %                             subplot(4,2,m)
                    %                             plot(para.bins,psthCT1(m).psth,'LineWidth', 2,'color','b');
                    %                             hold on;
                    %                             plot(para.bins,psthCT2(m).psth,'LineWidth', 2,'color','r');
                    %                             hold on;
                    %                             line([0 0],[0 max_resp], 'LineWidth', 1, 'Color', 'k')
                    %                             line([cue cue],[0 max_resp],'Color','k')
                    %                             line([cue+.5 cue+.5],[0 max_resp],'Color','k','linestyle','- -')
                    %                             line([sample sample],[0 max_resp],'Color','k')
                    %                             line([sample+.5 sample+.5],[0 max_resp],'Color','k','linestyle','- -')
                    %                             line([target target],[0 max_resp],'Color','k')
                    %                             axis([0 4.5 0 max_resp])
                    %
                    %                             if ~ isempty(psthCT2(m).psth)
                    %                                 subplot(4,2,m+4)
                    %                                 plot(para.bins,psthCT2(m).psth-psthCT1(m).psth,'LineWidth', 2,'color','k');
                    %                                 line([0 0],[-max_resp max_resp], 'LineWidth', 1, 'Color', 'k')
                    %                                 line([cue cue],[-max_resp max_resp],'Color','k')
                    %                                 line([cue+.5 cue+.5],[-max_resp max_resp],'Color','k','linestyle','- -')
                    %                                 line([sample sample],[-max_resp max_resp],'Color','k')
                    %                                 line([sample+.5 sample+.5],[-max_resp max_resp],'Color','k','linestyle','- -')
                    %                                 line([target target],[-max_resp max_resp],'Color','k')
                    %                                 axis([0 4.5 -max_resp max_resp ])
                    %                             end
                    %                         end
                end
                
            end
            %                 print('-djpeg', ['C:\work\Data_Analysis\APM_Data\DA_Data\stimulating\psthNeurons\SigIncreaseAny\' filename ])
            clf;
            
        end
        psthT1=[]; psthT2=[];
        if ~isempty(psthM)% & isempty(find(ntrStimT<2))
            disp(filename);
            for m = 1 : length(psthM)
                if normalize
                    %                     eval(['psthMall' num2str(m) '(n,:)=psthM(m).psth./maxpsth(n,1);'])
                    %                     eval(['psthNMall' num2str(m) '(n,:)=psthNM(m).psth./maxpsth(n,1);'])
                    if m<11
                        eval(['psthMall' num2str(m) '(n,:)=psthM(m).psth./fixrateMean(1);'])
                        eval(['psthNMall' num2str(m) '(n,:)=psthNM(m).psth./fixrateMean(1);'])
                    else
                        eval(['psthMall' num2str(m) '(n,:)=psthM(m).psth./fixrateMean(2);'])
                        eval(['psthNMall' num2str(m) '(n,:)=psthNM(m).psth./fixrateMean(2);'])
                    end
                else
                    if ~isempty(psthM(m).psth)
                        eval(['psthMall' num2str(m) '(n,:)=psthM(m).psth;'])
                        psthT1 = [psthT1; psthM(m).psth];
                    else
                        eval(['psthMall' num2str(m) '(n,:)=nan;'])
                    end
                    if ~isempty(psthNM(m).psth)
                        eval(['psthNMall' num2str(m) '(n,:)=psthNM(m).psth;'])
                        psthT2 = [psthT2; psthNM(m).psth];
                    else
                        eval(['psthNMall' num2str(m) '(n,:)=NaN;'])
                    end
                end
            end
            allpsth1(n,:)= nanmean(psthT1);%no stim
            allpsth2(n,:)= nanmean(psthT2);% stim
        else
            for m = 1 : 20
                eval(['psthMall' num2str(m) '(n,1:length(bin_edges))=NaN;'])
                eval(['psthNMall' num2str(m) '(n,1:length(bin_edges))=NaN;'])
            end
            %             allpsthC1(n,:)= psthCT1;
            %             allpsthC1(n,:)= psthCT2;
        end
        clear psthM psthNM;
    end
    %     end
end

% average fixrate through match/nonmatch if were seperate.
% for Nsheet = 1 : Npsth
%     eval(['fixrate' num2str(Nsheet)  '= nanmean(fixrate' num2str(Nsheet) '_1,2);']);
% end

%delete not enough classes
index =[];
for n = 1:20
    index =[ index find(isnan(eval(['psthMall' num2str(n) '(:,1)'])))'];
end

colors1={[0.678 0.922 1], [1 0.6 0.784]};
figure;hold on;
for n = 1 : 20
    lengends=[];
    nn=0;
    if n ==6 | n==11 | n==16
        figure;hold on;
    end
    
    for m = n
        nn = nn+1;
        
        %         NormalizeMax = max((eval(['psthMall' num2str(n)])),[],2);
        %         NormalizeMax=NormalizeMax';
        %         eval(['psthM(m,:)= nanmean(psthMall' num2str(n) ')./NormalizeMax;']);
        %         NormalizeMax = max((eval(['psthNMall' num2str(n)])),[],2);
        %         eval(['psthNM(m,:)= nanmean(psthNMall' num2str(n) ')./NormalizeMax;'])
        eval(['psthM(m,:)= nanmean(psthMall' num2str(n) ');']);
        eval(['psthNM(m,:)= nanmean(psthNMall' num2str(n) ');'])
        eval(['[ indexNostim a] = find(psthMall' num2str(n) '(:,1)>=0);']);
        eval(['[ indexStim a] = find(psthNMall' num2str(n) '(:,1)>=0);']);
        sizeNostim = length(indexNostim);
        sizeStim = length(indexStim);
        eval(['s = std(psthMall' num2str(n) '(indexNostim,:))./sqrt(sizeNostim);']);
        eval(['s1 = std(psthNMall'  num2str(n) '(indexStim,:))./sqrt(sizeStim);']);
        
        if mod(n,5)>0
            subplot(5,1,mod(n,5))
        else
            subplot(5,1,5)
        end
        hold on
        plot(para.bins,psthM(m,:),colors{1},'linewidth',2);
        plot(para.bins,psthNM(m,:),colors{2},'linewidth',2);
        lengends{1}=['noStim-' classDef{n}];
        lengends{2}=['Stim-' classDef{n}];
        legend(lengends);
        shadedplot(para.bins,psthM(m,:)-s, psthM(m,:)+s,colors1{1},colors1{1});
        %         alpha(0.5);
        shadedplot(para.bins,psthNM(m,:)-s1, psthNM(m,:)+s1,colors1{2},colors1{2});
        %         alpha(0.5);
        
        plot(para.bins,psthM(m,:),colors{1},'linewidth',2);
        plot(para.bins,psthNM(m,:),colors{2},'linewidth',2);
        %        if isSeperateMNM
        %            eval(['psthNM(m,:) = nanmean(psthNMall' num2str(m) '1);'])
        %            plot(para.bins,psthNM(m,:),colors{nn},'linewidth',2);
        %        end
    end    
    line([cue cue],[0 100],'Color','b')
    line([cue+.5 cue+.5],[0 100],'Color','b')
    line([sample sample],[0 100],'Color','b')
    line([sample+.5 sample+.5],[0 100],'Color','b')
    line([target target],[0 100],'Color','b')
    axis(para.axisNumber)
    title(figname);
    %     hold off
end


%  plot all psth together
allpsthMean1= nanmean(allpsth1);
allpsthMean2= nanmean(allpsth2);
nostimNtr =length( find(~isnan(allpsth1(:,1))));
stimNtr =length( find(~isnan(allpsth2(:,1))));
psthData =[];
para.legends =  {['NonStim =' num2str(nostimNtr)],['Stim =' num2str(stimNtr)]};
% for Nsheet = 1 : Npsth
eval(['psthData =[allpsthMean1; allpsthMean2;' '];']);

%     eval(['para.fixrate =[para.fixrate; nanmean(fixrate' num2str(Nsheet) ',1)];']);
% end
hb1 = figure;
hold on;
% bins = bin_edges+0.5*bin_width;
for n = 1 : size(psthData,1)
    eval(['plot(para.bins,psthData(n,:),colors{n} , ''LineWidth'', 2);']);
end
legend(para.legends);
hold on;
line([cue cue],[0 100],'Color','k','LineStyle','--')
line([cue+.5 cue+.5],[0 100],'Color','k','LineStyle','--')
line([sample sample],[0 100],'Color','k','LineStyle','--')
line([sample+.5 sample+.5],[0 100],'Color','k','LineStyle','--')
line([target target],[0 100],'Color','k','LineStyle','--')
axis(para.axisNumber)
xlabel('time (s)','fontsize',10,'fontweight','b');
title(['PSTH nonStim vs stim '])
axis on



% cue plot

for Nsheet = 1 : Npsth
    eval(['psthM1= nanmean(psthCM1);']);
    eval(['psthM2 = nanmean(psthCM2);']);
    eval(['psthM3 = nanmean(psthCM3);']);
    eval(['psthM4 = nanmean(psthCM4);']);
    eval(['psthNM1= nanmean(psthCNM1);']);
    eval(['psthNM2 = nanmean(psthCNM2);']);
    eval(['psthNM3 = nanmean(psthCNM3);']);
    eval(['psthNM4 = nanmean(psthCNM4);']);
    eval(['psthM1= nanmean(psthCM1);']);
    eval(['psthM2 = nanmean(psthCM2);']);
    eval(['psthM3 = nanmean(psthCM3);']);
    eval(['psthM4 = nanmean(psthCM4);']);
    eval(['psthNM1= nanmean(psthCNM1);']);
    eval(['psthNM2 = nanmean(psthCNM2);']);
    eval(['psthNM3 = nanmean(psthCNM3);']);
    eval(['psthNM4 = nanmean(psthCNM4);']);
end

%White task: Max Classes:
psthData =[];
para.legends =  {'1stNonStim','1stStim'};
% for Nsheet = 1 : Npsth
eval(['psthData =[psthM1; psthNM1;' '];']);
%     eval(['para.fixrate =[para.fixrate; nanmean(fixrate' num2str(Nsheet) ',1)];']);
% end
hb1 = figure;
hold on;
% bins = bin_edges+0.5*bin_width;
for n = 1 : size(psthData,1)
    eval(['plot(para.bins,psthData(n,:),colors{n} , ''LineWidth'', 2);']);
end
legend(para.legends);
hold on;
line([cue cue],[0 100],'Color','k','LineStyle','--')
line([cue+.5 cue+.5],[0 100],'Color','k','LineStyle','--')
line([sample sample],[0 100],'Color','k','LineStyle','--')
line([sample+.5 sample+.5],[0 100],'Color','k','LineStyle','--')
line([target target],[0 100],'Color','k','LineStyle','--')
axis(para.axisNumber)
xlabel('time (s)','fontsize',10,'fontweight','b');
title([figname ' white [1:5]'])
axis on

%White task: opp Classes:
psthData =[];
para.legends =  {'1stNonStim','1stStim'};
% for Nsheet = 1 : Npsth
eval(['psthData =[psthM2; psthNM2;' '];']);
%     eval(['para.fixrate =[para.fixrate; nanmean(fixrate' num2str(Nsheet) ',1)];']);
% end
hb1 = figure;
hold on;
% bins = bin_edges+0.5*bin_width;
for n = 1 : size(psthData,1)
    eval(['plot(para.bins,psthData(n,:),colors{n} , ''LineWidth'', 2);']);
end
legend(para.legends);
hold on;
line([cue cue],[0 100],'Color','k','LineStyle','--')
line([cue+.5 cue+.5],[0 100],'Color','k','LineStyle','--')
line([sample sample],[0 100],'Color','k','LineStyle','--')
line([sample+.5 sample+.5],[0 100],'Color','k','LineStyle','--')
line([target target],[0 100],'Color','k','LineStyle','--')
axis(para.axisNumber)
xlabel('time (s)','fontsize',10,'fontweight','b');
title([figname ' white [6:10]'])
axis on


%White task: Max Classes:
psthData =[];
para.legends =  {'2ndNonStim','2ndStim'};
% for Nsheet = 1 : Npsth
eval(['psthData =[psthM3; psthNM3;' '];']);
%     eval(['para.fixrate =[para.fixrate; nanmean(fixrate' num2str(Nsheet) ',1)];']);
% end
hb1 = figure;
hold on;
% bins = bin_edges+0.5*bin_width;
for n = 1 : size(psthData,1)
    eval(['plot(para.bins,psthData(n,:),colors{n} , ''LineWidth'', 2);']);
end
legend(para.legends);
hold on;
line([cue cue],[0 100],'Color','k','LineStyle','--')
line([cue+.5 cue+.5],[0 100],'Color','k','LineStyle','--')
line([sample sample],[0 100],'Color','k','LineStyle','--')
line([sample+.5 sample+.5],[0 100],'Color','k','LineStyle','--')
line([target target],[0 100],'Color','k','LineStyle','--')
axis(para.axisNumber)
xlabel('time (s)','fontsize',10,'fontweight','b');
title([figname ' blue max [11:14]'])
axis on

%Blue task: opp Classes:
psthData =[];
para.legends =  {'2ndNonStim','2ndStim'};
% for Nsheet = 1 : Npsth
eval(['psthData =[psthM4; psthNM4;' '];']);
%     eval(['para.fixrate =[para.fixrate; nanmean(fixrate' num2str(Nsheet) ',1)];']);
% end
hb1 = figure;
hold on;
% bins = bin_edges+0.5*bin_width;
for n = 1 : size(psthData,1)
    eval(['plot(para.bins,psthData(n,:),colors{n} , ''LineWidth'', 2);']);
end
legend(para.legends);
hold on;
line([cue cue],[0 100],'Color','k','LineStyle','--')
line([cue+.5 cue+.5],[0 100],'Color','k','LineStyle','--')
line([sample sample],[0 100],'Color','k','LineStyle','--')
line([sample+.5 sample+.5],[0 100],'Color','k','LineStyle','--')
line([target target],[0 100],'Color','k','LineStyle','--')
axis(para.axisNumber)
xlabel('time (s)','fontsize',10,'fontweight','b');
title([figname ' blue opp [16:19]'])
axis on


% Align 2nd stimuli In RF
% white task:[4 6];Blue[14 16];
figure;hold on;
lengends=[];
psthM_SecInRFW=nanmean([psthM(4,:);psthM(6,:)]);%no Stim
psthNM_SecInRFW=nanmean([psthNM(4,:);psthNM(6,:)]);%Stim
hold on
plot(para.bins,psthM_SecInRFW,colors{1},'linewidth',2);
plot(para.bins,psthNM_SecInRFW,colors{2},'linewidth',2);

lengends{1}=['noStim-4 6'];
lengends{2}=['Stim-4 6'];
legend(lengends);
line([cue cue],[0 100],'Color','b')
line([cue+.5 cue+.5],[0 100],'Color','b')
line([sample sample],[0 100],'Color','b')
line([sample+.5 sample+.5],[0 100],'Color','b')
line([target target],[0 100],'Color','b')
axis(para.axisNumber)
title(figname);

% Align 2nd stimuli In RF
%Blue[14 16];
figure;hold on;
lengends=[];
psthM_SecInRFB=nanmean([psthM(14,:); psthM(15,:);psthM(16,:)]);%no Stim
psthNM_SecInRFB=nanmean([psthNM(14,:);psthNM(15,:);psthNM(16,:)]);%Stim
hold on
plot(para.bins,psthM_SecInRFB,colors{1},'linewidth',2);
plot(para.bins,psthNM_SecInRFB,colors{2},'linewidth',2);
lengends{1}=['noStim-14 15 16'];
lengends{2}=['Stim-14 15 16'];
legend(lengends);
line([cue cue],[0 100],'Color','b')
line([cue+.5 cue+.5],[0 100],'Color','b')
line([sample sample],[0 100],'Color','b')
line([sample+.5 sample+.5],[0 100],'Color','b')
line([target target],[0 100],'Color','b')
axis(para.axisNumber)
title(figname);


% Align 2nd stimuli out/opp RF
% white task:[19];Blue[11 19];
figure;hold on;
lengends=[];
psthM_SecInRFW=nanmean([psthM(1,:);psthM(9,:)]);%no Stim
psthNM_SecInRFW=nanmean([psthNM(1,:);psthNM(9,:)]);%Stim
hold on
plot(para.bins,psthM_SecInRFW,colors{1},'linewidth',2);
plot(para.bins,psthNM_SecInRFW,colors{2},'linewidth',2);

lengends{1}=['noStim-1 9'];
lengends{2}=['Stim-1 9'];
legend(lengends);
line([cue cue],[0 100],'Color','b')
line([cue+.5 cue+.5],[0 100],'Color','b')
line([sample sample],[0 100],'Color','b')
line([sample+.5 sample+.5],[0 100],'Color','b')
line([target target],[0 100],'Color','b')
axis(para.axisNumber)
title(figname);

% Align 2nd stimuli out/opp RF
%Blue[11 19 20];
figure;hold on;
lengends=[];
psthM_SecInRFB=nanmean([psthM(11,:);psthM(19,:);psthM(20,:)]);%no Stim
psthNM_SecInRFB=nanmean([psthNM(11,:);psthNM(19,:);psthNM(20,:)]);%Stim
hold on
plot(para.bins,psthM_SecInRFB,colors{1},'linewidth',2);
plot(para.bins,psthNM_SecInRFB,colors{2},'linewidth',2);
lengends{1}=['noStim-11 19 20'];
lengends{2}=['Stim-11 19 20'];
legend(lengends);
line([cue cue],[0 100],'Color','b')
line([cue+.5 cue+.5],[0 100],'Color','b')
line([sample sample],[0 100],'Color','b')
line([sample+.5 sample+.5],[0 100],'Color','b')
line([target target],[0 100],'Color','b')
axis(para.axisNumber)
title(figname);

% Align 2nd stimuli 90 out RF
% white task:[2  7];Blue[12 17];
figure;hold on;
lengends=[];
psthM_SecInRFW=nanmean([psthM(2,:);psthM(7,:)]);%no Stim
psthNM_SecInRFW=nanmean([psthNM(2,:);psthNM(7,:)]);%Stim
hold on
plot(para.bins,psthM_SecInRFW,colors{1},'linewidth',2);
plot(para.bins,psthNM_SecInRFW,colors{2},'linewidth',2);

lengends{1}=['noStim-2 7'];
lengends{2}=['Stim-2 7'];
legend(lengends);
line([cue cue],[0 100],'Color','b')
line([cue+.5 cue+.5],[0 100],'Color','b')
line([sample sample],[0 100],'Color','b')
line([sample+.5 sample+.5],[0 100],'Color','b')
line([target target],[0 100],'Color','b')
axis(para.axisNumber)
title(figname);

% Align 2nd stimuli In RF
%Blue[12 17];
figure;hold on;
lengends=[];
psthM_SecInRFB=nanmean([psthM(12,:); psthM(17,:)]);%no Stim
psthNM_SecInRFB=nanmean([psthNM(12,:);psthNM(17,:)]);%Stim
hold on
plot(para.bins,psthM_SecInRFB,colors{1},'linewidth',2);
plot(para.bins,psthNM_SecInRFB,colors{2},'linewidth',2);
lengends{1}=['noStim-12 17'];
lengends{2}=['Stim-12 17'];
legend(lengends);
line([cue cue],[0 100],'Color','b')
line([cue+.5 cue+.5],[0 100],'Color','b')
line([sample sample],[0 100],'Color','b')
line([sample+.5 sample+.5],[0 100],'Color','b')
line([target target],[0 100],'Color','b')
axis(para.axisNumber)
title(figname);



 %plot Null together
%[5 10 15 20];
figure;hold on;
lengends=[];
hold on
plot(para.bins,nanmean(psthM([5 10],:)),colors{1},'linewidth',2);
plot(para.bins,nanmean(psthNM([5 10],:)),colors{2},'linewidth',2);
plot(para.bins,nanmean(psthM([15 20],:)),colors{3},'linewidth',2);
plot(para.bins,nanmean(psthNM([15 20],:)),colors{4},'linewidth',2);

lengends{1}=['noStim-null-1st'];
lengends{2}=['Stim-null-first'];
lengends{3}=['noStim-null-2nd'];
lengends{4}=['Stim-null-2nd'];
legend(lengends);
line([cue cue],[0 100],'Color','b')
line([cue+.5 cue+.5],[0 100],'Color','b')
line([sample sample],[0 100],'Color','b')
line([sample+.5 sample+.5],[0 100],'Color','b')
line([target target],[0 100],'Color','b')
axis(para.axisNumber)
title(figname);


%
% %Blue task
% psthData =[];
% para.legends =  {'2ndNoStim','2ndStim'};
% for Nsheet = 1 : Npsth
%
%     eval(['psthData =[psthData; psth' num2str(Nsheet) '3;' '];']);
% %     eval(['para.fixrate =[para.fixrate; nanmean(fixrate' num2str(Nsheet) ',1)];']);
% end
% hb1 = figure;
% hold on;
% % bins = bin_edges+0.5*bin_width;
% for n = 1 : size(psthData,1)
%     eval(['plot(para.bins,psthData(n,:),colors{n} , ''LineWidth'', 2);']);
% end
% legend(para.legends);
% shadedplot(para.bins, psthData(1,:)-s11, psthData(1,:)+s12, colors1{1}, colors1{1});
% hold on
% shadedplot(para.bins, psthData(2,:)-s21, psthData(2,:)+s22, colors1{2}, colors1{2});
% for n = 1 : size(psthData,1)
%     hold on;
%     eval(['plot(para.bins,psthData(n,:),colors{n} , ''LineWidth'', 2);']);
% end
% hold on;
% line([cue cue],[0 100],'Color','k','LineStyle','--')
% line([cue+.5 cue+.5],[0 100],'Color','k','LineStyle','--')
% line([sample sample],[0 100],'Color','k','LineStyle','--')
% line([sample+.5 sample+.5],[0 100],'Color','k','LineStyle','--')
% line([target target],[0 100],'Color','k','LineStyle','--')
% axis(para.axisNumber)
% xlabel('time (s)','fontsize',10,'fontweight','b');
% title(['PSTH sigCue pfc nonSti vs stim Blue '])
% axis on
%
% % plot white blue together
% psthData =[];
% para.legends =  {'1st nonStim', '2nd nonStim','1st Stim','2nd Stim'};
% for Nsheet = 1 : Npsth
%       eval(['psthData =[psthData; psth' num2str(Nsheet) '1;' '];']);
%     eval(['psthData =[psthData; psth' num2str(Nsheet) '3;' '];']);
% %     eval(['para.fixrate =[para.fixrate; nanmean(fixrate' num2str(Nsheet) ',1)];']);
% end
% hb1 = figure;
% hold on;
% % bins = bin_edges+0.5*bin_width;
% for n = 1 : size(psthData,1)
%     eval(['plot(para.bins,psthData(n,:),colors{n} , ''LineWidth'', 2);']);
% end
% legend(para.legends);
% % shadedplot(para.bins, psthData(1,:)-s11, psthData(1,:)+s12, colors1{1}, colors1{1});
% % hold on
% % shadedplot(para.bins, psthData(2,:)-s21, psthData(2,:)+s22, colors1{2}, colors1{2});
% % for n = 1 : size(psthData,1)
% %     hold on;
% %     eval(['plot(para.bins,psthData(n,:),colors{n} , ''LineWidth'', 2);']);
% % end
% hold on;
% line([cue cue],[0 100],'Color','k','LineStyle','--')
% line([cue+.5 cue+.5],[0 100],'Color','k','LineStyle','--')
% line([sample sample],[0 100],'Color','k','LineStyle','--')
% line([sample+.5 sample+.5],[0 100],'Color','k','LineStyle','--')
% line([target target],[0 100],'Color','k','LineStyle','--')
% axis(para.axisNumber)
% xlabel('time (s)','fontsize',10,'fontweight','b');
% title(['PSTH sigCue pfc nostim vs stim Blue '])
% axis on