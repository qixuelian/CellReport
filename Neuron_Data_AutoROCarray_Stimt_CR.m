% This script automatically runs Neuron_Data_ROCarray_bootstrap using the list of
% significant neurons that is also used for PSTH plotting (after truncating,
% ANOVA2, and dropping color significant neurons).
% Change the values of TI, TW, and iter on Neuron_Data_ROCarray_bootstrap if necessary.
% All ROC values from each neuron and ROC low/high values from bootstrap
% analysis will be separately saved as .dat file.
% Import these dat files into separate excel sheets and average ROC values for each bin
% to perform population analysis.
% Last edited: Oct 05 2009
% By Fumi Katsuki
% Jan 12 2010, XQ
% define  TW and TI in MainScripts
% save data as .mat file
% save plot into .emf, .ai, .fig

clear;
warning off MATLAB:divideByZero
%%% Change Parameter %%%
Area = 'PFC3mon sig cue spatial'
filerange = 'Spatial'
IsSeperate = 1;
CueDelay = 0; % caculate according to Cue or Cue Delay firing rate
if CueDelay
    filerange = [ filerange 'CueDelay']
    phase ='cuedelay';
else
     filerange = [ filerange 'Cue']
    phase ='cuerate';
end
%%% change over here
AlignPhase= 4; % 1 cue; 2 -cudelay, 3 sample, 4 sampledelay
legend_txt = {'nonstim','stim'};
excelName = 'C:\work\DataBase\StimulationFilename_neuron_ODRdistVar.xlsx'%
sheetName1 = 'neuronNoSigDecreaseEnoughTrial' %'NoSigAnyChange'%'sigIncreaseAnyNoMixEnough'%'neuronNoSigDecreaseEnoughTrial' %'allNeurons' %'neuronNoSigDecreaseEnoughTrial'%'neuronNoSigchangeEnoughTrial'%'neuronSigIncreaseEnoughTrial'%'neuronNoSigDecreaseEnoughTrial'%'neuronUseSigIncreaseFix_use'%'neuronUseNoSigChangeFix'%'neuronUseIncreaseFix'%'neuronUseIncreaseFix'%'SigAnySigIncreaFiringrate' %'neurons_ROC' %%'mdlPFC_Elv';'pre_PFC3monSNRmax_spatial'
sheetName2 = 'post_PFC_3monkey'%%'mdlPFC_elv2';'ElvNeuronLeftPFCSpatial''post_PFC3monSNRmax_spatial''ScopElvSpatialdlPFC'

cWindow = [-1 5.5]; % XQ
% avoid timewindow protruding into other period, resulting in an artificial
% significant ROA
%Align spikes by CueOn, SampleOn,Target/Reward_On seperately
cWindow_temp =[-1 0;0 2;0 2;0 1.5];
TI = 0.05; % sliding bin size for ROC 0.05
TW = 0.1;  % integration time window to collect firing rates
iter = 100;
cutoff = 96; % p<0.05 5th percentile, for iter=100, cutoff=96, for iter=1000, cutoff=951
TWms = TW*1000
steps = TI* 1000
Bootstrap = iter;


% out put Data
% outPutDir = ['C:\Matlab2006B\work\Data_Analysis\APM_Data\DA_Data\ROC\'];
outPutDir =  ['C:\work\Data_Analysis\APM_Data\DA_Data\ROC\'];
outPutData =  [outPutDir 'ROC_',Area,'_NonSeperate',filerange,'TW',num2str(TWms),'TI', num2str(steps)];
outPutPlot =  [outPutDir 'ROCpop_',Area,'_NonSeperate',filerange,'TW',num2str(TWms), 'TI', num2str(steps)];

%%%%%%%%%%%%%%%%%%%%%%%%
%PreTraining
[Neurons_num Neurons_txt] = xlsread(excelName,sheetName1);
% [Neurons_num Neurons_txt] = xlsread('C:\Matlab2006B\work\DataBase\NeuronData_Hemisphere.xlsx',sheet1);
Neurons = [Neurons_txt(:,1) num2cell(Neurons_num(:,1))];
cuetime = cWindow(1):TI:cWindow(2); % xq
results_w(1:length(Neurons),1:length(cuetime))=NaN;
results_w_stim(1:length(Neurons),1:length(cuetime))=NaN;
results_b(1:length(Neurons),1:length(cuetime))=NaN;
results_b_stim(1:length(Neurons),1:length(cuetime))=NaN;
AUCss_w(1:length(Neurons),1:length(cuetime))=NaN;
AUCss_w_stim(1:length(Neurons),1:length(cuetime))=NaN;
AUCss_b(1:length(Neurons),1:length(cuetime))=NaN;
AUCss_b_stim(1:length(Neurons),1:length(cuetime))=NaN;

switch AlignPhase
    case 1 %cue
        cueStart = length(find(cuetime<=0));
        cueEnd = length(find(cuetime<=0.3));
        cuePoint = cueStart:cueEnd;
        plotEnd = length(find(cuetime<=3));
    case 2 % cuedelay
        
        cueStart = length(find(cuetime<=0.5));
        cueEnd = length(find(cuetime<=1.5));
        cuePoint = cueStart:cueEnd;
        plotEnd = length(find(cuetime<=3));
    case 3% sample
        cueStart = length(find(cuetime<=1.5));
        cueEnd = length(find(cuetime<=1.8));
        cuePoint = cueStart:cueEnd;
        plotEnd = length(find(cuetime<=3));
    case 4%sampledelay
        cueStart = length(find(cuetime<=2));
        cueEnd = length(find(cuetime<=3));
        cuePoint = cueStart:cueEnd;
        plotEnd = length(find(cuetime<=3));
end

invert = 1;
% fh1=figure;
for n = 1:length(Neurons)
    filename = [Neurons{n,1},'_',num2str(Neurons{n,2})];
    try
        %     [ROCarea_M ROCarea_NM ROCarea_M_NM ROCarea_NM_M] = Neuron_Data_ROCarray_bootstrap_XQ(filename,IsFeature,cWindow_temp,cWindow,TI,TW,iter,cutoff);
        if IsSeperate% filename,cWindow,TI,TW,iter,cutoff,phase,invert
            if AlignPhase >2
            [ROCarea_M ROCarea_NM ROCarea_M_NM ROCarea_NM_M] = Neuron_Data_ROCarray_bootstrap_odrdistVar_Stim_align2nd(filename,cWindow,TI,TW,iter,cutoff,phase,invert);
            % [ROCarea_M ROCarea_NM ROCarea_M_NM ROCarea_NM_M] = Neuron_Data_ROCarray_bootstrap_odrdistVar_Stim_align2ndTrial(filename,cWindow,TI,TW,iter,cutoff,phase,invert);
            else
                [ROCarea_M ROCarea_NM ROCarea_M_NM ROCarea_NM_M]  = Neuron_Data_ROCarray_bootstrap_odrdistVar_Stim(filename,cWindow,TI,TW,iter,cutoff,phase,invert);
            % [ROCarea_M ROCarea_NM ROCarea_M_NM ROCarea_NM_M]  = Neuron_Data_ROCarray_bootstrap_odrdistVar_Stim_HUA(filename,cWindow,TI,TW,iter,cutoff,phase,invert);
            % [ROCarea_M ROCarea_NM ROCarea_M_NM ROCarea_NM_M] = Neuron_Data_ROCarray_bootstrap_odrdistVar_Stim_throughTrial(filename,cWindow,TI,TW,iter,cutoff,phase,invert);
            end
            results_w(n,1:length(cuetime)) = ROCarea_M;%white
            results_w_stim(n,1:length(cuetime)) = ROCarea_NM;%white stim
            results_b(n,1:length(cuetime)) = ROCarea_M_NM;%blue
            results_b_stim(n,1:length(cuetime)) = ROCarea_NM_M;%blue stim
            %smoothed
%             AUCs1(n,1:length(cuetime))=smooth(ROCarea_M,3)';
            AUCss_w(n,1:length(cuetime))=smooth(ROCarea_M,5)';
%             AUCs1_stim(n,1:length(cuetime))=smooth(ROCarea_NM,3)';
            AUCss_w_stim(n,1:length(cuetime))=smooth(ROCarea_NM,5)';
%             AUCs2(n,1:length(cuetime))=smooth(ROCarea_M_NM,3)';
            AUCss_b(n,1:length(cuetime))=smooth(ROCarea_M_NM,5)';
%             AUCs2_stim(n,1:length(cuetime))=smooth(ROCarea_NM_M,3)';
            AUCss_b_stim(n,1:length(cuetime))=smooth(ROCarea_NM_M,5)';
            %             figure(fh1)
            %             clf;
            %             hold on;
            %             plot(cuetime, ROCarea_NM,'-b','LineWidth', 2);
            %             % legend('pre','post')
            %             line([0 0],[0 1],'Color','k')
            %             line([.5 .5],[0 1],'Color','k','LineStyle','--')
%             line([1.5 1.5],[0 1],'Color','k')
%             line([2 2],[0 1],'Color','k','LineStyle','--')
%             line([3 3],[0 1],'Color','k','LineStyle','--')
%             line([-2 7.5],[0.5 0.5],'color','k','linestyle','--')
%             axis([cWindow(1) cWindow(2) 0 1])
%             title([filename])
%             set(fh1,'Color',[.5 .5 .5])
%             saveas(gcf,['C:\work\Data_Analysis\APM_Data\DA_Data\stimulating\ROC\filePlot\ROCarea_NM_' filename],'jpeg')
        else
            [ROCarea_M] = Neuron_Data_ROCarray_bootstrap_XQ_Nonseperate(filename,IsFeature,cWindow,TI,TW,iter,cutoff,CueDelay);
            results_M(n,1:length(cuetime)) = ROCarea_M;
            AUCs(n,1:length(cuetime))=smooth(ROCarea_M,3)';
            AUCss(n,1:length(cuetime))=smooth(ROCarea_M,5)';
        end
    catch
        if IsSeperate
            results_w(n,1:length(cuetime)) =nan;
            results_w_stim(n,1:length(cuetime)) = nan;
             results_b(n,1:length(cuetime)) =nan;
            results_b_stim(n,1:length(cuetime)) = nan;
             AUCss_w(n,1:length(cuetime))=nan;
            AUCss_w_stim(n,1:length(cuetime))=nan;
             AUCss_b(n,1:length(cuetime))=nan;
            AUCss_b_stim(n,1:length(cuetime))=nan; 
%             AUCs3(n,1:length(cuetime))=nan;
%             AUCss3(n,1:length(cuetime))=nan;
%              AUCs4(n,1:length(cuetime))=nan;
%             AUCss4(n,1:length(cuetime))=nan;
        else
            results_M(n,1:length(cuetime)) =nan;
            AUCs(n,1:length(cuetime))=nan;
            AUCss(n,1:length(cuetime))=nan;
        end
        disp([lasterr '  on:  ' filename]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TWms = TW*1000
steps = TI* 1000
Bootstrap = iter;

%%% Population Analysis %%%
% Plot nonstim 
popresults = nanmean(results_w);
figure;
hold on;
plot(cuetime,popresults,'-b','LineWidth', 2);
% plot(cuetime,nanmean(results_M2),'-g','LineWidth', 2);
plot(cuetime,nanmean(results_w_stim),'-r','LineWidth', 2);
% legend('pre','post')
legend(legend_txt);
line([0 0],[0 1],'Color','k')
line([.5 .5],[0 1],'Color','k','LineStyle','--')
line([1.5 1.5],[0 1],'Color','k')
line([2 2],[0 1],'Color','k','LineStyle','--')
line([3 3],[0 1],'Color','k','LineStyle','--')
line([-2 7.5],[0.5 0.5],'color','k','linestyle','--')
axis([cWindow(1) cWindow(2) 0.2 0.8])
title(['ROCpop','  ',Area,' ',filerange,'white (TW:',num2str(TWms),'ms'])
% print('-dmeta', [outPutPlot 'MatchMatch.emf'])
% print('-dill', [outPutPlot 'MatchMatch.ai'])
% saveas(gcf,[outPutPlot 'MatchMatch.fig'])
% Plot NonMatch vs NonMatch
popresults = nanmean(results_b);
figure;
hold on;
plot(cuetime,popresults,'-b','LineWidth', 2);
% plot(cuetime,nanmean(results_NM2),'-g','LineWidth', 2);
plot(cuetime,nanmean(results_b_stim),'-r','LineWidth', 2);
% legend('pre','post')
legend(legend_txt);
line([0 0],[0 1],'Color','k')
line([.5 .5],[0 1],'Color','k','LineStyle','--')
line([1.5 1.5],[0 1],'Color','k')
line([2 2],[0 1],'Color','k','LineStyle','--')
line([3 3],[0 1],'Color','k','LineStyle','--')
line([-2 7.5],[0.5 0.5],'color','k','linestyle','--')
axis([cWindow(1) cWindow(2) 0.2 0.8])
title(['ROCpop','  ',Area,' ',filerange,'blue (TW:',num2str(TWms),'ms'])
% print('-dmeta', [outPutPlot 'NM_NM.emf'])
% print('-dill', [outPutPlot 'NM_NM.ai'])
% saveas(gcf,[outPutPlot 'NM_NM.fig'])

%%
clims = [0.5 1];
meanCue=nanmean(results_w(:,cuePoint ),2); %% cue 21:30  31:51
index=find(~isnan(meanCue));
[SequentialPosition,SequentialTrialIndex]=sortrows(meanCue(index));
SequentialTrialIndex=index(SequentialTrialIndex);
% [SequentialPosition,SequentialTrialIndex]=sortrows(Position,-1);% '-'sort Poke interval descend %SequentialTrialIndex
for m=1:size(SequentialTrialIndex,1)%size(results_w,1)
    results_w_plot(m,:)=AUCss_w(SequentialTrialIndex(m),1:plotEnd);
%     results_w_plot(m,:)=results_w(SequentialTrialIndex(m),1:plotEnd);
end

% plot the figure
figure('color',[1 1 1])
imagesc(results_w_plot,clims)
xlabel('Time from fixation (s)','fontsize',12)
ylabel('Number of neuron','fontsize',12)
% set(gca, 'XTick', [0.5,10:10:80]);
 set(gca, 'XTick', [1,25:25:250]);
set(gca,'XTickLabel',{'-1','-0.5','0','0.5','1','1.5','2','2.5','3'})
colorbar('ylim',[0.5 1], 'YTick', [0.5:0.1:1],'position',[0.92 0.11 0.03 0.815],'fontsize',10)
title(['white task no stim: '  ])

%%
meanCue=nanmean(results_w_stim(:,cuePoint),2); %% cue 21:30
index=find(~isnan(meanCue));
[SequentialPosition,SequentialTrialIndex]=sortrows(meanCue(index));
SequentialTrialIndex=index(SequentialTrialIndex);
% [SequentialPosition,SequentialTrialIndex]=sortrows(Position,-1);% '-'sort Poke interval descend %SequentialTrialIndex
for m=1:size(SequentialTrialIndex,1)%size(results_w_stim,1)
    results_w_stim_plot(m,:)=AUCss_w_stim(SequentialTrialIndex(m),1:plotEnd);
%     results_w_stim_plot(m,:)=results_w_stim(SequentialTrialIndex(m),1:plotEnd);
end

% plot the figure
figure('color',[1 1 1])
imagesc(results_w_stim_plot,clims)
xlabel('Time from fixation (s)','fontsize',12)
ylabel('Number of neuron','fontsize',12)
% set(gca, 'XTick', [0.5,10:10:80]);
 set(gca, 'XTick', [1,25:25:250]);
set(gca,'XTickLabel',{'-1','-0.5','0','0.5','1','1.5','2','2.5','3'})
colorbar('ylim',[0.5 1], 'YTick', [0.5:0.1:1],'position',[0.92 0.11 0.03 0.815],'fontsize',10)
title(['white task stim: '  ])

%%
meanCue=nanmean(results_b(:,cuePoint),2); %% cue 21:30
index=find(~isnan(meanCue));
[SequentialPosition,SequentialTrialIndex]=sortrows(meanCue(index));
SequentialTrialIndex=index(SequentialTrialIndex);

% [SequentialPosition,SequentialTrialIndex]=sortrows(Position,-1);% '-'sort Poke interval descend %SequentialTrialIndex
for m=1:size(SequentialTrialIndex,1)
    results_b_plot(m,:)=AUCss_b(SequentialTrialIndex(m),1:plotEnd);
%      results_b_plot(m,:)=results_b(SequentialTrialIndex(m),1:plotEnd);
end

% plot the figure
figure('color',[1 1 1])
imagesc(results_b_plot,clims)
xlabel('Time from fixation (s)','fontsize',12)
ylabel('Number of neuron','fontsize',12)
% set(gca, 'XTick', [0.5,10:10:80]);
 set(gca, 'XTick', [1,25:25:250]);
set(gca,'XTickLabel',{'-1','-0.5','0','0.5','1','1.5','2','2.5','3'})
colorbar('ylim',[0.5 1], 'YTick', [0.5:0.1:1],'position',[0.92 0.11 0.03 0.815],'fontsize',10)
title(['blue task no stim: '  ])

%%
meanCue=nanmean(results_b_stim(:,cuePoint),2); %% cue 21:30 cd 31:51
index=find(~isnan(meanCue));
[SequentialPosition,SequentialTrialIndex]=sortrows(meanCue(index));
SequentialTrialIndex=index(SequentialTrialIndex);
% [SequentialPosition,SequentialTrialIndex]=sortrows(Position,-1);% '-'sort Poke interval descend %SequentialTrialIndex
for m=1:size(SequentialTrialIndex,1)
    results_b_stim_plot(m,:)=AUCss_b_stim(SequentialTrialIndex(m),1:plotEnd);
%     results_b_stim_plot(m,:)=results_b_stim(SequentialTrialIndex(m),1:plotEnd);
end

% plot the figure
figure('color',[1 1 1])
imagesc(results_b_stim_plot,clims)
xlabel('Time from fixation (s)','fontsize',12)
ylabel('Number of neuron','fontsize',12)
% set(gca, 'XTick', [0.5,10:10:80]);
 set(gca, 'XTick', [1,25:25:250]);
set(gca,'XTickLabel',{'-1','-0.5','0','0.5','1','1.5','2','2.5','3'})
colorbar('ylim',[0.5 1], 'YTick', [0.5:0.1:1],'position',[0.92 0.11 0.03 0.815],'fontsize',10)
title(['blue task stim: '  ])