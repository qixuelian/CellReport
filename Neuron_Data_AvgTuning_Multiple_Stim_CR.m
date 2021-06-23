% function betas = Neuron_Data_Tuning(varargin)

% 2 inputs
% betas = Neuron_Data_Tuning('filename_neuron#',auto_var)
% This function inputs a compressed spatial data file and an automating
% switch and outputs either a 3D figure or the gaussian values of the
% spatial receptive field 'betas' which equals: [baseline, amplitude,
% X, Y, standard deviation].  The X and Y coordinates represent the center
% of the receptive field, standard deviation is the width.  Passing
% auto_var = 1 automates the processes and only outputs the gaussian
% numbers, auto_var = 0 displays the figure.
% 1 input
% Same thing but it just displays the figure
% 5-16-06
% Mar-08-2010
% XQ
%OCT-14-2010 XQ
% creat from Neuron_Data_AvgTuning_XQ
%for spatial task only; align max firing locations, not hemispheric

clear
auto_var = 1;
phase ='cue' %'cue'% 'cuedelay';  % cue or cuedelay
Ndata = 1;
excell_filename = 'C:\work\DataBase\StimulationFilename_neuron_ODRdistVar.xlsx';
sheetName1 = 'neuronNosigDecreaseEnoughTrial' % 'neuronUseSigdecreaseFix' %'neuronNosigDecreaseEnoughTrial' %neuronUseSigdecreaseFix'sigIncreaseAnyNoMixEnough'%'neuronNoSigDecreaseEnoughTrial'

% monkey = {'elv' 'adr' 'ben'};monkeys = ['elv' 'adr' 'ben'];
colorsNum = {[0 0 1]  [1 0 0]};
colors = {'b' 'r'};
legend_txt={'pre' 'post'};
task = 'spatial';
Area = ['dlPFC'];
FigTitle = [sheetName1 ' ' phase];
reverseN = [6 7 8 9 10 1 2 3 4 5 16 17 18 19 20 11 12 13 14 15];
cue_dur =0.5;
cue_dur2=0.2;%use rate 
cuedelay_dur =1

for Nsheet = 1
    [Neurons_num Neurons_txt] = xlsread(excell_filename,eval(['sheetName' num2str(Nsheet)]));
    warning off MATLAB:divideByZero
    Neurons = [Neurons_txt(:,1) num2cell(Neurons_num(:,1))];
    class_means = ones(9,size(Neurons,1)) * NaN;
    class_meansG =ones(10,size(Neurons,1)) * NaN;
    ntr = 1;
    for n = 1:length(Neurons)
        filename = [Neurons{n,1},'_',num2str(Neurons{n,2})];
        load(filename)
        if ~isempty(MatData)
            varCue = [];varCueStim=[];  nIndex=[];
            fix = [];
            for j = 1:length(MatData.class)
                for nnn = 1 : length(MatData.class(j).ntr)
                    TS_2 = MatData.class(j).ntr(nnn).TS;
%                     MatData.class(j).ntr(nnn).cuerate_200 =length(find((TS_2>MatData.class(j).ntr(nnn).Cue_onT) & (TS_2<=(MatData.class(j).ntr(nnn).Cue_onT + cue_dur))))/cue_dur;
                    MatData.class(j).ntr(nnn).samplerate_200 =length(find((TS_2>MatData.class(j).ntr(nnn).Cue_onT+ cue_dur+cuedelay_dur) & (TS_2<=(MatData.class(j).ntr(nnn).Cue_onT + cue_dur+cuedelay_dur+cue_dur2))))/cue_dur2;
                end
                
                if strcmpi( phase,'cuedelay')
                    cuerate(j) = mean([MatData.class(j).ntr.cuedelay]);
                    nIndex = find([MatData.class(j).ntr.Stim] == 1);
                    varCueStim(j) = nanmean([MatData.class(j).ntr(nIndex).sampledelay]);
                    nIndex = find([MatData.class(j).ntr.Stim] == 0);
                    varCue(j)  = nanmean([MatData.class(j).ntr(nIndex).sampledelay]);
                else
                    cuerate(j) = mean([MatData.class(j).ntr.cuerate]);
                    nIndex = find([MatData.class(j).ntr.Stim] == 1);
                    varCueStim(j) = nanmean( [MatData.class(j).ntr(nIndex).samplerate_200]);
                    nIndex = find([MatData.class(j).ntr.Stim] == 0);
                    varCue(j) = nanmean( [MatData.class(j).ntr(nIndex).samplerate_200]);
                end
                if isfield(MatData.class(j).ntr,'fixrate')
                    fix(j) = mean([MatData.class(j).ntr.fixrate]);
                else
                    fix(j) = mean([MatData.class(j).ntr.fix]);
                end
            end
            fixrate = nanmean(fix);
            
            [a max_class] = max([nanmean(cuerate(1:5)) nanmean(cuerate(6:10)) nanmean(cuerate(11:14)) nanmean(cuerate(16:19))]);
            if mod(max_class,2)==0;
                isReverse ='r';
                varCueStim =varCueStim(reverseN);
                varCue=varCue(reverseN);
            else
                isReverse ='n';
            end
            CueStimW(5,ntr)=nanmean(varCueStim([4 6]));%max
            CueStimW(4,ntr)=varCueStim(3);%45 from max
            CueStimW(3,ntr)=nanmean(varCueStim([2 7]));%90
            CueStimW(2,ntr)=varCueStim(8);%135
            CueStimW(1,ntr)=nanmean(varCueStim([1 9]));%180 from max
            CueStimW(6,ntr)=varCueStim(3);%45 from max
            CueStimW(7,ntr)=nanmean(varCueStim([2 7]));%90
            CueStimW(8,ntr)=varCueStim(8);%135
            CueStimW(9,ntr)=nanmean(varCueStim([1 9]));%180 from max
            
            
            CueW(5,ntr)=nanmean(varCue([4 6]));%max
            CueW(4,ntr)=varCue(3);%45 from max
            CueW(3,ntr)=nanmean(varCue([2 7]));%90
            CueW(2,ntr)=varCue(8);%135
            CueW(1,ntr)=nanmean(varCue([1 9]));%180 from max
            CueW(6,ntr)=varCue(3);%45 from max
            CueW(7,ntr)=nanmean(varCue([2 7]));%90
            CueW(8,ntr)=varCue(8);%135
            CueW(9,ntr)=nanmean(varCue([1 9]));%180 from max
            
            CueStimB(5,ntr)=nanmean(varCueStim([14 16]));%max
            CueStimB(4,ntr)=varCueStim(13);%45 from max
            CueStimB(3,ntr)=nanmean(varCueStim([12 17]));%90
            CueStimB(2,ntr)=varCueStim(18);%135
            CueStimB(1,ntr)=nanmean(varCueStim([11 19]));%180 from max
            CueStimB(6,ntr)=varCueStim(13);%45 from max
            CueStimB(7,ntr)=nanmean(varCueStim([12 17]));%90
            CueStimB(8,ntr)=varCueStim(18);%135
            CueStimB(9,ntr)=nanmean(varCueStim([11 19]));%180 from max
            
            CueB(5,ntr)=nanmean(varCue([14 16]));%max
            CueB(4,ntr)=varCue(13);%45 from max
            CueB(3,ntr)=nanmean(varCue([12 17]));%90
            CueB(2,ntr)=varCue(18);%135
            CueB(1,ntr)=nanmean(varCue([11 19]));%180 from max
            CueB(6,ntr)=varCue(13);%45 from max
            CueB(7,ntr)=nanmean(varCue([12 17]));%90
            CueB(8,ntr)=varCue(18);%135
            CueB(9,ntr)=nanmean(varCue([11 19]));%180 from max
            
            ntr = ntr +1;
        else
            disp(['empty file ' filename]);
        end

    end

end
% end
% create matrix with average stimulus spikes, location, average fixation spikes
disp(excell_filename)

%white task
figure;
hold on;    xgrid=1:0.1:9;
for Nsheet = 1 : Ndata
%     [RW2 sW mW semW betasW modelFunW] =  gaus8loc_fit(eval(['class_meansWG'  num2str(Nsheet)]),colors);
%     [R2B sB mB semB betasB modelFunB] =  gaus8loc_fit(eval(['class_meansBG'  num2str(Nsheet)]),colors);
    [RW2 sW mW semW betasW modelFunW] =  gaus8loc_fit(CueStimW,colors);
    [R2B sB mB semB betasB modelFunB] =  gaus8loc_fit(CueW,colors);
    
    %Plot
    rate_plotW = modelFunW(betasW,xgrid);
%     errorbar(sW,mW,semW,'marker','o','linestyle','none','color',colorsNum{2})%stim
    plot(xgrid,rate_plotW,colors{2});
    hold on
    scatter([1:9], nanmean(CueStimW,2)','marker','o')
     scatter([1:9], nanmean(CueW,2)','marker','o')
    rate_plotB = modelFunB(betasB,xgrid);
%     errorbar(sB,mB,semB,'marker','o','linestyle','none','color',colorsNum{1})%no stim
    plot(xgrid,rate_plotB,colors{1});
end
xlim([0.5 10.5]);
%     legend(lengend_txt);
title([FigTitle 'white']);
%Blue task
figure;hold on;xgrid=1:0.1:9;
[RW2 sW mW semW betasW modelFunW] =  gaus8loc_fit(CueStimB,colors);
    [R2B sB mB semB betasB modelFunB] =  gaus8loc_fit(CueB,colors);
    
    %Plot
    rate_plotW = modelFunW(betasW,xgrid);
    errorbar(sW,mW,semW,'marker','o','linestyle','none','color',colorsNum{2})%(:,2)-m)
    plot(xgrid,rate_plotW,colors{2});
    hold on
    rate_plotB = modelFunB(betasB,xgrid);
    errorbar(sB,mB,semB,'marker','o','linestyle','none','color',colorsNum{1})%(:,2)-m)
    plot(xgrid,rate_plotB,colors{1});
xlim([0.5 10.5]);
%     legend(lengend_txt);
title([FigTitle ' blue']);

figure;hold on;xgrid=1:0.1:9;
[RW2 sW mW semW betasW modelFunW] =  gaus8loc_fit([ CueStimW CueStimB],colors);
    [R2B sB mB semB betasB modelFunB] =  gaus8loc_fit([CueW CueB],colors);
    
    %Plot
    rate_plotW = modelFunW(betasW,xgrid);
    errorbar(sW,mW,semW,'marker','o','linestyle','none','color',colorsNum{2})%(:,2)-m)
    plot(xgrid,rate_plotW,colors{2});
    hold on
    rate_plotB = modelFunB(betasB,xgrid);
    errorbar(sB,mB,semB,'marker','o','linestyle','none','color',colorsNum{1})%(:,2)-m)
    plot(xgrid,rate_plotB,colors{1});
xlim([0.5 10.5]);
%     legend(lengend_txt);
title([FigTitle ' white&blue']);

reg2=[ones(1,5) ;1:5]';

for n = 1 : size(CueStimW,2)
     aa=regress(CueW(1:5,n),reg2); 
     regressData(n,1)=aa(2);
    aa=regress(CueStimW(1:5,n),reg2); 
    regressData(n,2)=aa(2);
    aa=regress(CueB(1:5,n),reg2);
    regressData(n,3)=aa(2);
    aa=regress(CueStimB(1:5,n),reg2);
    regressData(n,4)=aa(2);
    
  end


% %SNR (m5-m1)/std1
% SNR1 = (nanmean(class_meansG1(5,:))-nanmean(class_meansG1(1,:)))./std(class_meansG1(1,:));
%     SNR2 = (nanmean(class_meansG2(5,:))-nanmean(class_meansG2(1,:)))./std(class_meansG2(1,:));

%line Plot
% [m,p] = grpstats(class_means1',{'1','2','3','4','5','6','7','8','9','10'},{'mean','sem'})
% n = length(m);
% figure;
% errorbar((1:n)',m,p,'linestyle','none','color',colors)%(:,2)-m)

% figure;
% rates = [real_means(6) real_means(7) real_means(8);real_means(5) real_means(9) real_means(1); real_means(4) real_means(3) real_means(2)];
% [X,Y]=meshgrid(-10:10:10, -10:10:10);
% [XI YI]=meshgrid(-10:1:10, -10:1:10);
% rateint=interp2(X,Y,rates,XI,YI);
% contourf(XI,YI,rateint,40)
% shading flat
