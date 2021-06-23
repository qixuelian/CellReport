function [ROCarea_1 ROCarea_1_stim ROCarea_2 ROCarea_2_stim] = Neuron_Data_ROCarray_bootstrap_odrdistVar_Stim_align2nd(filename,cWindow,TI,TW,iter,cutoff,phase,invert)

% Neuron_Data_ROCarray_bootstrap_XQ_odrdistVar_Stim% Time-resolved ROC analysis from Saliency experiment.
% Computes statistical significance of ROC values based on a bootstrap test
%
% Trials from arrays of different colors are pooled together.
% TI, TW parameters define size and overlap of sliding windows.
%
% Modified to make it compatible with 4 location tasks
% with some empty classes (e.g.NoJump4_D15_16+).
% By Fumi Katsuki & her advisor, Christos Constantinidis, Ph.D.
% Edited: Oct 03 2009
%
% Add lines to save figures at the end.
% By Fumi Katsuki
% Edited: Oct 05 2009
% align 2nd stimuli

oppClasslist = [5 6 7 8 1 2 3 4];% Spatial Class only XQ
load([filename '.mat']);
try
    MatData;
catch
    MatData = MUAData;
end
%%% Find the max response class %%%
try
    if ~isempty(MatData)
        %%%% find the class psth will be extracted from
        fixrate =[];
        for m = 1 :10;% length(MatData.class)
            %     cuerate(m) = nanmean([MatData.class(m).ntr.cuedelay]);
            fixrate = [fixrate [MatData.class(m).ntr.fix]];
        end
        fixMean(1,1) = nanmean(fixrate);
        
        fixrate =[];
        for m = 11 :20;% length(MatData.class)
            %     cuerate(m) = nanmean([MatData.class(m).ntr.cuedelay]);
            fixrate = [fixrate [MatData.class(m).ntr.fix]];
        end
        fixMean(1,2) = nanmean(fixrate);
        
        cuedelayT1=[];
        for m = [1:4]
            eval(['cuedelayT1 =[cuedelayT1 [MatData.class(m).ntr.' phase ']];'])
        end
        cuedelayT2=[];
        for m = [6:9]
            eval(['cuedelayT2 =[cuedelayT2 [MatData.class(m).ntr.'  phase ']];'])
        end
        
        cuedelayT3=[];
        for m = [11:14]
            eval(['cuedelayT3 =[cuedelayT3 [MatData.class(m).ntr.' phase ']];'])
        end
        cuedelayT4=[];
        for m = [16:19]
            eval(['cuedelayT4 =[cuedelayT4 [MatData.class(m).ntr.' phase ']];'])
        end
        
        [a b] = max([ nanmean(cuedelayT1)  nanmean(cuedelayT2)  nanmean(cuedelayT3)  nanmean(cuedelayT4)]);
        if mod(b,2)==0;
            isReverse ='r';
            if invert >0
                invert =1;
            else
                invert =0;
            end
        else
            isReverse ='n';
            invert =0;
        end
        if strcmpi(isReverse,'r')
            % if isReverse
            reverseN = [6 7 8 9 10 1 2 3 4 5 16 17 18 19 20 11 12 13 14 15];
        end
        
        if fixMean < 10000
            %             % White
            %             max_class1 = 1;% max-opp
            %             opp_class1 = 6;% opp-max;
            %             opp_class2 = 9;% max-max
            %             max_class2 = 4;%opp-opp
            %             % Blue
            %             max_class3 = 11;% max-opp
            %             opp_class3= 16;% opp-max;
            %             max_class4= 14;% max-max
            %             opp_class4 = 19;%opp-opp
            
            %%% Extract cue rate for the max_class_corr %%%
            ROCarea_1 = [];%white
            ROCarea_1_stim = [];
            ROCarea_2 = [];%blue
            ROCarea_2_stim = [];
            %             ROCarea_3 = [];
            %             ROCarea_4 = [];
            n=0;
            if  length(MatData.class) >= 20%max_class2 &  length([MatData.class(max_class2).ntr]) > 2
                for t = (cWindow(1)-TW/2)/TI:1:(cWindow(2)-TW/2)/TI
                    maxcl_cueall = [];
                    maxcl_cueall2 = [];
                    maxcl_cueall_stim = [];
                    maxcl_cueall2_stim = [];
                    n=n+1;
                    ntr_stim=0;
                    for m = [4 6]
                        max_class1 = m;
                        for i = 1:length(MatData.class(max_class1).ntr)
                            ntr_TS1 = MatData.class(max_class1).ntr(i).TS;
                            Cue_onT1 = MatData.class(max_class1).ntr(i).Cue_onT;
                            maxcl_TS1 = ntr_TS1(find(ntr_TS1 >= Cue_onT1+(TI*t) & ntr_TS1 <Cue_onT1+TI*t+TW));
                            if MatData.class(max_class1).ntr(i).Stim
                                maxcl_cueall_stim = [maxcl_cueall_stim length(maxcl_TS1)]; % Best Match
                                ntr_stim = ntr_stim+1;
                            else
                                maxcl_cueall = [maxcl_cueall length(maxcl_TS1)]; % Best Match
                            end
                        end
                    end
                    ntr_stim2 =0;
                    for m = [14 15 16]%11:14
                        max_class2=m;
                        for  i = 1:length(MatData.class(max_class2).ntr)
                            ntr_TS2 = MatData.class(max_class2).ntr(i).TS;
                            Cue_onT2 = MatData.class(max_class2).ntr(i).Cue_onT;
                            maxcl_TS2 = ntr_TS2(find(ntr_TS2 >= Cue_onT2+(TI*t) & ntr_TS2 < Cue_onT2+TI*t+TW));
                            if MatData.class(max_class2).ntr(i).Stim
                                maxcl_cueall2_stim = [maxcl_cueall2_stim length(maxcl_TS2)]; % Best NonMatch
                                ntr_stim2 = ntr_stim2+1;
                            else
                                maxcl_cueall2 = [maxcl_cueall2 length(maxcl_TS2)]; % Best NonMatch
                            end
                        end
                    end
                    %%% Extract cue rate for the opposite corr classes %%%
                    oppcl_cueall = [];
                    oppcl_cueall_stim = [];
                    oppcl_cueall2 = [];
                    oppcl_cueall2_stim = [];
                    ntr_stim3=0;
                    for m = [1 9] %6:9
                        opp_class1 = m;
                        for i = 1:length(MatData.class(opp_class1).ntr)
                            Cue_onT_opp1 = MatData.class(opp_class1).ntr(i).Cue_onT;
                            ntr_TS_opp1 = MatData.class(opp_class1).ntr(i).TS;
                            oppcl_TS1 = ntr_TS_opp1(find(ntr_TS_opp1 >= Cue_onT_opp1+(TI*t) & ntr_TS_opp1 < Cue_onT_opp1+TI*t+TW));
                            if MatData.class(opp_class1).ntr(i).Stim
                                oppcl_cueall_stim = [oppcl_cueall_stim length(oppcl_TS1)]; % OPP Match
                                ntr_stim3 = ntr_stim3+1;
                            else
                                oppcl_cueall = [oppcl_cueall length(oppcl_TS1)]; % OPP Match
                            end
                        end
                    end
                    ntr_stim4=0;
                    for m = [11 19 20]%16:19
                        opp_class2 = m;
                        for i = 1:length(MatData.class(opp_class2).ntr)
                            ntr_TS_opp2 = MatData.class(opp_class2).ntr(i).TS;
                            Cue_onT_opp2 = MatData.class(opp_class2).ntr(i).Cue_onT;
                            oppcl_TS2 = ntr_TS_opp2(find(ntr_TS_opp2 >= Cue_onT_opp2+(TI*t) & ntr_TS_opp2 < Cue_onT_opp2+TI*t+TW));
                            if MatData.class(opp_class2).ntr(i).Stim
                                ntr_stim4= ntr_stim4+1;
                                oppcl_cueall2_stim = [oppcl_cueall2_stim length(oppcl_TS2)]; % Opp NonMatch
                            else
                                oppcl_cueall2 = [oppcl_cueall2 length(oppcl_TS2)]; % Opp NonMatch
                            end
                        end
                    end
                    if ntr_stim >3 & ntr_stim2>3 & ntr_stim3>3 & ntr_stim4 >3
                        if ~strcmpi(isReverse,'r')
                            ROCarea_1(n) = arrayROC(maxcl_cueall,oppcl_cueall); % best match vs opp match
                            ROCarea_1_stim(n) = arrayROC(maxcl_cueall_stim,oppcl_cueall_stim); % best nonmatch vs opp NonMatch
                            ROCarea_2(n) = arrayROC(maxcl_cueall2,oppcl_cueall2); % best match vs opp match
                            ROCarea_2_stim(n) = arrayROC(maxcl_cueall2_stim,oppcl_cueall2_stim); % best nonmatch vs opp NonMatch
                        else
                            ROCarea_1(n) = arrayROC(oppcl_cueall,maxcl_cueall); % best match vs opp match
                            ROCarea_1_stim(n) = arrayROC(oppcl_cueall_stim,maxcl_cueall_stim); % best nonmatch vs opp NonMatch
                            ROCarea_2(n) = arrayROC(oppcl_cueall2,maxcl_cueall2); % best match vs opp match
                            ROCarea_2_stim(n) = arrayROC(oppcl_cueall2_stim,maxcl_cueall2_stim); % best nonmatch vs opp NonMatch
                            
                        end
                        alltrials=[maxcl_cueall oppcl_cueall];
                        totaltr_maxcl=length(maxcl_cueall);
                        totaltr_oppcl=length(oppcl_cueall);
                        ROCarea_boot=zeros(1,iter);
                        for it=1:iter
                            % Rearrange trials randomly
                            newindex=randperm(totaltr_maxcl+totaltr_oppcl);
                            maxcl_boot=zeros(1,totaltr_maxcl);
                            oppcl_boot=zeros(1,totaltr_oppcl);
                            for tr=1:totaltr_maxcl
                                maxcl_boot(tr)=alltrials(newindex(tr));
                            end
                            for tr=1:totaltr_oppcl
                                oppcl_boot(tr)=alltrials(newindex(tr+totaltr_maxcl));
                            end
                            % Compute ROC based on shuffled trials
                            ROCarea_boot1(it)=arrayROC(maxcl_boot,oppcl_boot);
                        end
                        % Compute Average and Standard Deviation of shuffled ROC value
                        %     ROCarea_est(n)=mean(ROCarea_boot);
                        %     ROCarea_std(n)=std(ROCarea_boot);
                        %     ROCarea_lo(n)=ROCarea_est(n) - 2*ROCarea_std(n);
                        %     ROCarea_hi(n)=ROCarea_est(n) + 2*ROCarea_std(n);
                        ROCarea_srt = sort(ROCarea_boot);
                        ROCarea_lo(n)=ROCarea_srt(iter-cutoff);
                        ROCarea_hi(n)=ROCarea_srt(cutoff);
                        
                    else
                        ROCarea_1 = nan;
                        ROCarea_1_stim=nan;
                        ROCarea_2 = nan;
                        ROCarea_2_stim=nan;
                        return
                    end
                end
            else
                ROCarea_1 = nan* ones( length((cWindow(1)-TW/2)/TI:1:(cWindow(2)-TW/2)/TI),1);
                ROCarea_1_stim=ROCarea_1;
                ROCarea_2 = ROCarea_1;
                ROCarea_2_stim=ROCarea_1;
                return
                %                 ROCarea_M = nan; % best match vs opp match
                %                 ROCarea_NM = nan; % best nonmatch vs opp NonMatch
                %                 ROCarea_M_NM(n) = nan; % best match vs opp NonMatch
                %                 ROCarea_NM_M(n) = nan; %best NonMatch vs opp match
            end
        else
            ROCarea_1 = nan;
            ROCarea_1_stim=nan;
            ROCarea_2 = nan;
            ROCarea_2_stim=nan;
            return
            %             ROCarea_M = nan; % best match vs opp match
            %             ROCarea_NM = nan; % best nonmatch vs opp NonMatch
            %             ROCarea_M_NM(n) = nan; % best match vs opp NonMatch
            %             ROCarea_NM_M(n) = nan; %best NonMatch vs opp match
        end
    end
catch
    lasterr
end

