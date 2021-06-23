% function Neuron_Data_PSTH_ODR_3intervals(Excel_neurons,sheetname)
% Plots Population PSTH for ODR task, providing detail for fixation period
% Identifies trials with at least 200 ms of pre-fixation period
% Also plots last second of intertrial interval
% 22-Oct-2016, Christos Constantinidis, PhD

%MIX White and blue tasks
%USE in paper

clear
Excel_neurons='C:\work\DataBase\StimulationFilename_neuron_ODRdistVar.xlsx' %'C:\work\DataBase\StimulationFilename_neuron_ODRdistVar';Stimulation_waveform
sheetname= 'neuronNoSigDecreaseEnoughTrial'
% 'allNeurons_use' % 'sigdecreaseAllNeuron' %'sigIncreaseAllNeuron' %'NochangeAllNeuron'  %'sigdecreaseAllNeuron' %'noSigNeurons_Use'% 'neuronsigDecreseenoughtrial' %'neuronNoSigDecreaseEnoughTrial'
% 'neuronNoSigDecreaseEnoughgru' 'neuronUseNoSigChangeFix'%sigIncreaseAnyNoMixEnough'%'neuronUseSigdecreaseFix'%'neuronUseNoSigChangeFix'%'neuronUseSigIncreaseFix_use'%'neuronNoSigchangeEnoughTrial'
[Neurons_num Neurons_txt] = xlsread(Excel_neurons,sheetname);
warning off MATLAB:divideByZero
Neurons = [Neurons_txt(:,1) num2cell(Neurons_num(:,1))];

timeend = 2.5;
bin_width =0.025;% 0.025;  % 50 milliseconds bin
bin_edges=-1:bin_width:timeend;   % -0.8
bins = bin_edges+0.5*bin_width;
m_counter=0;
allTS=[];
fixTS=[];
interTS=[];
BehTemp=['hec184_2','hec189_2','hec190_2','hec193_2','hec194_2','hec195_2','hec196_2'];
disp(length(Neurons))
stim =1;
for n = 1:length(Neurons)
    Profilename = [Neurons{n,1},'_',num2str(Neurons{n,2})];
    Behfilename = lower([Neurons{n,1}]);
    filename = lower([Neurons{n,1}]);
    m_counter=0;    allTS=[];    fixTS=[];    interTS=[];
    m_counter_stim=0;    allTS_stim=[];    fixTS_stim=[];  interTS_stim=[];
    try
        load(Profilename); % load spike mat datafile
        load(Behfilename); % load behavioral datafile
        if strfind(BehTemp,Behfilename)
            load([Behfilename(1:7) '1CATends']);
            if length(task_ends)>=3
                disp(['multiple Beh']);
                %         MatDataT = MatData;
                tt = AllData;
                for nn= 3 : length(task_ends)
                    %             ttype=nn;
                    try
                        load ([filename(1:7) num2str(nn)]);
                        if ~isfield( tt.trials,'StimOnTime')
                            tt.trials(1).StimOnTime=[];
                        end
                        if ~isfield( AllData.trials,'StimOnTime')
                            AllData.trials(1).StimOnTime=[];
                        end
                        tt.trials= orderfields(tt.trials,AllData.trials);
                        
                        tt.trials(task_ends(nn-1)-task_ends(1)+1:task_ends(nn-1)-task_ends(1)+length(AllData.trials))=AllData.trials;
                        %             .trialnum +task_ends(ttype)-task_ends(ttype-1)+1;
                        %matdata
                        for nclass = 1 : length(MatData.class)
                            diffNumtrial=diff([MatData.class(nclass).ntr.trialnum]);
                            indexdiff = find(diffNumtrial<=0)+1;
                            if ~isempty(indexdiff)
                                if length(indexdiff)>1
                                    endofNtr = indexdiff(2);
                                else
                                    endofNtr = length(MatData.class(nclass).ntr);
                                end
                                for nnn = indexdiff(1) : endofNtr%length(MatData.class(nclass).ntr)
                                    origTrialNum = MatData.class(nclass).ntr(nnn).trialnum -1+task_ends(nn-1)-task_ends(nn);
                                    MatData.class(nclass).ntr(nnn).trialnum  = origTrialNum + task_ends(nn-1)-task_ends(1);
                                end
                            end
                        end
                    end
                end
                AllData =tt;
            end
        end
        
        %white/remeber-first task and blue task together
        ntr = 0;rate_inter=[]; rate_inter_stim =[];ntr_stim=0;
        rate_fix=[]; rate_fix_stim =[];
        for class_num = 1:length(MatData.class)
            stim=0;
            index =find([MatData.class(class_num).ntr(:).Stim] == stim);
             MatData1.class(class_num).ntr = [];
            MatData1.class(class_num).ntr = MatData.class(class_num).ntr(index);
            for m = 1:length(MatData1.class(class_num).ntr)
                try
                    tn = MatData1.class(class_num).ntr(m).trialnum;
                    % real fixation onset time
                    Fix_onT = AllData.trials(tn).FixOn - AllData.trials(tn).time;
                    % real end of trial - 200 ms
                    End_onT = length(AllData.trials(tn).EyeData)/500 - 0.2;
                    m_counter = m_counter + 1;
                    ntr = ntr +1;
                    TS = MatData1.class(class_num).ntr(m).TS - MatData1.class(class_num).ntr(m).Cue_onT;
                    allTS = [allTS TS];
                    TS = MatData1.class(class_num).ntr(m).TS - Fix_onT;
                    fixTS = [fixTS TS];
                    rate_fix(ntr)= length(find(TS>0 & TS<=0.5));
                    TS = MatData1.class(class_num).ntr(m).TS - End_onT;
                    interTS=[interTS TS];                    
                    rate_inter(ntr)= length(find(TS>-1 & TS<=0));
                catch
                    disp([lasterr ' on ' Profilename])
                end
            end
            %with stimulation
            stim=1;
            index =find([MatData.class(class_num).ntr(:).Stim] == stim);
            MatData2.class(class_num).ntr =[];
            MatData2.class(class_num).ntr = MatData.class(class_num).ntr(index);
            %             rate_inter_stim =[];ntr_stim=0;
            for m = 1:length(MatData2.class(class_num).ntr)
                try
                    tn = MatData2.class(class_num).ntr(m).trialnum;
                    % real fixation onset time
                    Fix_onT = AllData.trials(tn).FixOn - AllData.trials(tn).time;
                    % real end of trial - 200 ms
                    End_onT = length(AllData.trials(tn).EyeData)/500 - 0.2;
                     m_counter_stim = m_counter_stim + 1;
                    ntr_stim=ntr_stim+1;
                    TS = MatData2.class(class_num).ntr(m).TS - MatData2.class(class_num).ntr(m).Cue_onT;
                    allTS_stim = [allTS_stim TS];
                    TS = MatData2.class(class_num).ntr(m).TS - Fix_onT;
                    fixTS_stim = [fixTS_stim TS];
                    rate_fix_stim(ntr_stim)= length(find(TS>0 & TS<=0.5));
                    TS = MatData2.class(class_num).ntr(m).TS - End_onT;
                    interTS_stim=[interTS_stim TS];                   
                    rate_inter_stim(ntr_stim)= length(find(TS>-1 & TS<=0));
                     
                catch
                    disp([lasterr ' on ' Profilename])
                end
            end
        end
        if m_counter_stim>=1
            ntrs = m_counter;
            Propsth_w(n,1:length(bin_edges)) =histc(allTS,bin_edges)/(bin_width*ntrs);
            Fixpsth_w(n,1:length(bin_edges)) =histc(fixTS,bin_edges)/(bin_width*ntrs);
            Interpsth_w(n,1:length(bin_edges)) =histc(interTS,bin_edges)/(bin_width*ntrs);
            
            ntrs = m_counter_stim;
            Propsth_stim_w(n,1:length(bin_edges)) =histc(allTS_stim,bin_edges)/(bin_width*ntrs);
            Fixpsth_stim_w(n,1:length(bin_edges)) =histc(fixTS_stim,bin_edges)/(bin_width*ntrs);
            Interpsth_stim_w(n,1:length(bin_edges)) =histc(interTS_stim,bin_edges)/(bin_width*ntrs);
            rate_inter_mean(n,:) = [nanmean(rate_inter) nanmean(rate_inter_stim)];
            rate_fix_mean(n,:) = [nanmean(rate_fix) nanmean(rate_fix_stim)];
            [h p(n,1)]= ttest2(rate_inter_stim,rate_inter);
            [h p_fix(n,1)]= ttest2(rate_fix_stim,rate_fix);
        else
            Propsth_w(n,1:length(bin_edges)) =nan;
            Fixpsth_w(n,1:length(bin_edges)) =nan;
            Interpsth_w(n,1:length(bin_edges)) =nan;
            Propsth_stim_w(n,1:length(bin_edges)) =nan;
            Fixpsth_stim_w(n,1:length(bin_edges)) =nan;
            Interpsth_stim_w(n,1:length(bin_edges)) = nan;
            rate_inter_mean(n,1:2) = nan;rate_fix_mean(n,1:2) = nan;
            p(n,1)= nan; p_fix(n,1)= nan;
        end
        
    catch
        disp([lasterr 'on file ' Profilename])
    end
end
definepsthmax=60;
figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot white/remeber-first task
Fixpsth=Fixpsth_w;
Interpsth=Interpsth_w;
Propsth=Propsth_w;
Fixpsth_stim=Fixpsth_stim_w;
Interpsth_stim=Interpsth_stim_w;
Propsth_stim=Propsth_stim_w;

colors1 = {[0.678 0.922 1], [1 0.6 0.784]};

subplot(3,2,3)
FixpsthM = nanmean(Fixpsth);
FixpsthM_stim = nanmean(Fixpsth_stim);

[ indexNostim a] = find(Fixpsth(:,1)>=0);
[ indexStim a] = find(Fixpsth_stim(:,1)>=0);
sizeNostim = length(indexNostim);
sizeStim = length(indexStim);
s = std(Fixpsth(indexNostim,:))./sqrt(sizeNostim);
s1 = std(Fixpsth_stim(indexStim,:))./sqrt(sizeStim);

% for k=1:length(bin_edges)
%     PsthSEM(k) = std(Fixpsth(:,k))/sqrt(n);
% end
% PropsthH=FixpsthM + PsthSEM;
% PropsthL=FixpsthM - PsthSEM;
% shadedplot(bins,PropsthH, PropsthL,[0.7 0.7 0.7],[0.7 0.7 0.7]);
hold on
shadedplot(bins,FixpsthM-s, FixpsthM+s,colors1{1},colors1{1});
alpha(0.5);
shadedplot(bins,FixpsthM_stim-s1, FixpsthM_stim+s1,colors1{2},colors1{2});
alpha(0.5);
plot(bins,FixpsthM,'b','LineWidth',1)%,'linestyle','- -');
plot(bins,FixpsthM_stim,'r','LineWidth',1)%,'linestyle','- -');
axis([0 .5 5 definepsthmax])
box('off')
% xlabel('Time s')
% ylabel('Firing Rate spikes/s')

subplot(3,2,5)
InterpsthM = nanmean(Interpsth);
InterpsthM_stim = nanmean(Interpsth_stim);

[ indexNostim a] = find(Interpsth(:,1)>=0);
[ indexStim a] = find(Interpsth_stim(:,1)>=0);
sizeNostim = length(indexNostim);
sizeStim = length(indexStim);
s = std(Interpsth(indexNostim,:))./sqrt(sizeNostim);
s1 = std(Interpsth_stim(indexStim,:))./sqrt(sizeStim);

% for k=1:length(bin_edges)
%     PsthSEM(k) = std(Interpsth(:,k))/sqrt(n);
% end
% PropsthH=InterpsthM + PsthSEM;
% PropsthL=InterpsthM - PsthSEM;
% shadedplot(bins,PropsthH, PropsthL,[0.7 0.7 0.7],[0.7 0.7 0.7]);
hold on
shadedplot(bins,InterpsthM-s, InterpsthM+s,colors1{1},colors1{1});
alpha(0.5);
shadedplot(bins,InterpsthM_stim-s1, InterpsthM_stim+s1,colors1{2},colors1{2});
alpha(0.5);
plot(bins,InterpsthM,'b','LineWidth',1)%,'linestyle','- -');
plot(bins,InterpsthM_stim,'r','LineWidth',1)%,'linestyle','- -');
line([0 0], [0 60],'color','k')
axis([-1 0 5 definepsthmax])
% xlabel('Time s')
% ylabel('Firing Rate spikes/s')
box('off')

subplot(3,1,1)
PropsthM = nanmean(Propsth);
PropsthM_stim = nanmean(Propsth_stim);

[ indexNostim a] = find(Propsth(:,1)>=0);
[ indexStim a] = find(Propsth_stim(:,1)>=0);
sizeNostim = length(indexNostim);
sizeStim = length(indexStim);
s = std(Propsth(indexNostim,:))./sqrt(sizeNostim);
s1 = std(Propsth_stim(indexStim,:))./sqrt(sizeStim);
% for k=1:length(bin_edges)
%     PsthSEM(k) = std(Propsth(:,k))/sqrt(n);
% end
% PropsthH=PropsthM + PsthSEM;
% PropsthL=PropsthM - PsthSEM;
%
% shadedplot(bins,PropsthH, PropsthL,[0.7 0.7 0.7],[0.7 0.7 0.7]);
hold on
shadedplot(bins,PropsthM-s, PropsthM+s,colors1{1},colors1{1});
alpha(0.5);
shadedplot(bins,PropsthM_stim-s1, PropsthM_stim+s1,colors1{2},colors1{2});
alpha(0.5);
plot(bins,PropsthM,'b','LineWidth',1)%,'linestyle','- -');
plot(bins,PropsthM_stim,'r','LineWidth',1)%,'linestyle','- -');
line([0 0], [0 60],'color','k')
line([2 2], [0 60],'color','k')
axis([-0.5 2.5 5 definepsthmax])
% xlabel('Time s')
% ylabel('Firing Rate spikes/s')
box('off')
title([sheetname '  w + b'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %plot blue/remeber-second task
% figure; hold on;
% Fixpsth=Fixpsth_b;
% Interpsth=Interpsth_b;
% Propsth=Propsth_b;
%
% Fixpsth_stim=Fixpsth_stim_b;
% Interpsth_stim=Interpsth_stim_b;
% Propsth_stim=Propsth_stim_b;
%
% subplot(3,2,3)
% FixpsthM = nanmean(Fixpsth);
% FixpsthM_stim = nanmean(Fixpsth_stim);
%
% % for k=1:length(bin_edges)
% %     PsthSEM(k) = std(Fixpsth(:,k))/sqrt(n);
% % end
% % PropsthH=FixpsthM + PsthSEM;
% % PropsthL=FixpsthM - PsthSEM;
% % shadedplot(bins,PropsthH, PropsthL,[0.7 0.7 0.7],[0.7 0.7 0.7]);
% hold on
% plot(bins,FixpsthM,'b','LineWidth',1)%,'linestyle','- -');
% plot(bins,FixpsthM_stim,'r','LineWidth',1)%,'linestyle','- -');
%
% axis([0 .5 5 definepsthmax])
% box('off')
% % xlabel('Time s')
% % ylabel('Firing Rate spikes/s')
%
% subplot(3,2,5)
% InterpsthM = nanmean(Interpsth);
% InterpsthM_stim = nanmean(Interpsth_stim);
% % for k=1:length(bin_edges)
% %     PsthSEM(k) = std(Interpsth(:,k))/sqrt(n);
% % end
% % PropsthH=InterpsthM + PsthSEM;
% % PropsthL=InterpsthM - PsthSEM;
% % shadedplot(bins,PropsthH, PropsthL,[0.7 0.7 0.7],[0.7 0.7 0.7]);
% hold on
% plot(bins,InterpsthM,'b','LineWidth',1)%,'linestyle','- -');
% plot(bins,InterpsthM_stim,'r','LineWidth',1)%,'linestyle','- -');
% line([0 0], [0 60],'color','k')
% axis([-1 0 5 definepsthmax])
% % xlabel('Time s')
% % ylabel('Firing Rate spikes/s')
% box('off')
%
% subplot(3,1,1)
% PropsthM = nanmean(Propsth);
% PropsthM_stim = nanmean(Propsth_stim);
% % for k=1:length(bin_edges)
% %     PsthSEM(k) = std(Propsth(:,k))/sqrt(n);
% % end
% % PropsthH=PropsthM + PsthSEM;
% % PropsthL=PropsthM - PsthSEM;
% %
% % shadedplot(bins,PropsthH, PropsthL,[0.7 0.7 0.7],[0.7 0.7 0.7]);
% hold on
% plot(bins,PropsthM,'b','LineWidth',1)%,'linestyle','- -');
% plot(bins,PropsthM_stim,'r','LineWidth',1)%,'linestyle','- -');
% line([0 0], [0 60],'color','k')
% line([2 2], [0 60],'color','k')
% axis([-0.5 2.5 5 definepsthmax])
% % xlabel('Time s')
% % ylabel('Firing Rate spikes/s')
% box('off')
% title([sheetname '  blue'])



