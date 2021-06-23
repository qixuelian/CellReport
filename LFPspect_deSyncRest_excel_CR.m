% function LFPspect_1(fileName)
%This function creates power spectra for LFP data as obtained from the APM
%module in Constantinidis' lab.  It uses a high-pass filter. It calls
%plot_LFP_nofilt.m
%no photodiode signals

%read excelto decide which rials to use and when the start point

clear
%parameters:
winTotal = 180;      %record time in second
winC = 4;       % step window
winStart = 1 % start time point-- not in use, use in excel settings

%we use 1000 ms to avoid spill-over effect from cue presentation
winMS = 500;      %match stimulus time window in ms
fpass = [0 100];  %frequency band to be displayed in spectrum ([fmin fmax])
%grulfp_1_CH2 grudesync1_ch2 'hecdesync_12_ch2' 'hecdesync_13_ch2' 'hecdesync_14_ch2'
[Neurons_num Neurons_txt] = xlsread('C:\work\DataBase\stimulationfilename_desync.xlsx','gruuse')%Sheet1'hecuse');%Neurons_ken.xlsx Neurons_DDRGrayVar.xlsx','Neurons_ODR');%Neurons_Grub.xls','ODR';Neurons_Lemon.xls','temp'

% fileName = {'hecdesync_14_ch2'}%{'hecdesync_5_ch2' 'hecdesync_6_ch2' 'hecdesync_7_ch2' 'hecdesync_9_ch2' ...
    %'hecdesync_10_ch2' 'hecdesync_11_ch2'}%'hecdesync_11_ch2'%'hecdesync_11_ch2' %'hecdesync_9_ch2'%;'hecdesync_6_ch2' grudesync1_ch2
    %small {'hecdesync_6_ch2' 'hecdesync_10_ch2'  }
 trials_excludes_list = {[]}%{[17] [] [17] [] [7 15 16] [2:3] };
 
%   LFP_Data_Plot_deSync(fileName{1})
  

%parameters for fft
pad =0;
% tapers = [2 3];

 %high-pass filter w/cutoff 4 Hz:
%hifilt=fir1(250, 2*4/500, 'high', hamming(251));  %linear-phase method,
%uses signal processing toolbox NOT GOOD.
% [Bfilt, Afilt] = butter(5, 2*4/500, 'high'); %IIR filter (non-linear phase) BETTER!
[Bfilt, Afilt] = butter(5, 2*0.001/500, 'high'); %IIR filter (non-linear phase) BETTER!
catchlist = [];

% fileName(end - 3:end) = '';% for hecdesync9_ch2.mat
% if ~iscell(fileName)
    load('hecdesync_14_ch2');
%     trials_use = [2: 14 16 18:length(LFPData.class(1).ntr)-1];% [2: 14 16:length(LFPData.class(1).ntr)];%[8: length(LFPData.class(1).ntr)-1];
% else
%    
%     for n =  1: length(fileName)
%         file = fileName{n};
%         load(file);
%         t = 1 :length(LFPData.class(1).ntr);
% %         trials_use = find(~ismember(trials_use,trials_excludes_list{n});
%         %         t.class(1).ntr()=LFPData.class(1).ntr(trials_use);
%     end
%     
% end

SP = [];
CU = [];
% DE = [];
% MS = [];

% try
%     Accepted = LFP_Trial_Check(fileName);
TotalAccepted = 30;
%%% for first trial collecting each 6 seconds
recordTime = 180; % seconds
SampleRate = LFPData.class(1).ntr(1).Sampling_Rate;
Con = ceil(recordTime*LFPData.class(1).ntr(1).Sampling_Rate);            %data point for 6 seconds
baseline_1 =[];
TempLFP = LFPData.class(1).ntr(1).LFP(1:Con); % very first trial
nb = recordTime/6;
TempLFP = reshape(TempLFP,6*SampleRate ,nb)';

% TempLFP = TempLFP';
% nb = length(TempLFP);
for k = 1:size(TempLFP, 1)  %convolve raw trace with filter
    TempLFP(k, :) = filtfilt(Bfilt, Afilt, TempLFP(k, :));
end
BaseLFP=TempLFP;
totaltrial=[];
%

for n = 1: size(Neurons_txt,1)
    clear LFPData
    file = [Neurons_txt{n,1} '_ch' num2str(Neurons_num(n,1))];
    disp(sprintf('Processing %s ...', file));
    load(file);
    t = 2 :length(LFPData.class(1).ntr);
     trials_use = eval(Neurons_txt{n,3});
     winStart=Neurons_num(n,3);
     totaltrial=[totaltrial trials_use];
%     trials_use =t((~ismember(t,trials_excludes_list{n})));
    if TotalAccepted > 20% 30 pre para = 30
        for i = 1:length(LFPData.class)
            for j = trials_use %
                %                 if Accepted.class(i).ntr(j)
                Con = ceil(6*SampleRate);            %data point for 6 seconds
                Clast = ceil(180*SampleRate)-1;
                S = LFPData.class(i).ntr(j).LFP(winStart.*SampleRate+1 :winStart.*SampleRate+Con);  %voltage trace during first 6 seconds (mV)
                C = LFPData.class(i).ntr(j).LFP(Clast-Con+1:Clast);  %voltage trace during  last 6 seconds(mV)
                
                S = S';
                nb = length(S);
                for k = 1:size(S, 1)  %convolve raw trace with filter
                    S(k, :) = filtfilt(Bfilt, Afilt, S(k, :));
                end
                if ~isempty(SP)
                    [a, b] = size(SP);
                    if b < nb,
                        SP = [SP NaN*ones(a, nb - b)];
                    elseif b > nb,
                        S = [S NaN*ones(1,b - nb)];
                    end
                end
                SP = [SP;S];
                
                
                C = C';
                nb = length(C);
                for k = 1:size(C, 1)  %convolve raw trace with filter
                    C(k, :) = filtfilt(Bfilt, Afilt, C(k, :));
                end
                if ~isempty(CU)
                    [a, b] = size(CU);
                    if b < nb,
                        CU = [CU NaN*ones(a, nb - b)];
                    elseif b > nb,
                        C = [C NaN*ones(1,b - nb)];
                    end
                end
                CU = [CU;C];
                
                %                 end
            end
        end
    end
end
[totalpowerF,powF,frS,varS,ntrials]=getpower3_2(SP,SampleRate);
[totalpowerC,powC,frS,varC,ntrials]=getpower3_2(CU,SampleRate);
[totalpowerB,powB,frS,varB,ntrials]=getpower3_2(BaseLFP,SampleRate);


%plot from fr [1. 20]
% USE in paper
a =  find(frS<=1);
b=  find(frS<=20);
figure;
% hold on
% semilogy(frS, totalpowerB, 'g-', 'linewidth', 2)
% legend('first 6 sec', 'last 6 sec');
semilogy(frS(a(end)+1:5:b(end)+1), totalpowerF(a(end)+1:5:b(end)+1), 'b-', 'linewidth', 2)
hold on
semilogy(frS(a(end)+1:5:b(end)+1), totalpowerC(a(end)+1:5:b(end)+1), 'r-', 'linewidth', 2)
legend('first 6 sec', 'last 6 sec')
ha = shadedplot_semilogy(frS(a(end)+1:5:b(end)+1),totalpowerF(a(end)+1:5:b(end)+1).*10.^(-varS(a(end)+1:5:b(end)+1)),  totalpowerF(a(end)+1:5:b(end)+1).*10.^varS(a(end)+1:5:b(end)+1),[1 0.7 0.7],[1 0.7 0.7] ); %first area is red
hold on
hb = shadedplot_semilogy(frS(a(end)+1:5:b(end)+1),  totalpowerC(a(end)+1:5:b(end)+1).*10.^(-varC(a(end)+1:5:b(end)+1)),  totalpowerC(a(end)+1:5:b(end)+1).*10.^(varC(a(end)+1:5:b(end)+1)),[0.7 0.7 1],[0.7 0.7 1]); 
hold on
semilogy(frS(a(end)+1:5:b(end)+1), totalpowerF(a(end)+1:5:b(end)+1), 'r-', 'linewidth', 2)
hold on
semilogy(frS(a(end)+1:5:b(end)+1), totalpowerC(a(end)+1:5:b(end)+1), 'b-', 'linewidth', 2)
ax = axis;
axis([2 20 ax(3) ax(4)])
xlim([1 20]);
ylim([50000 3000000]);
set(gca, 'XTick', fpass(1):10:fpass(2));
xlabel('Frequency (Hz)')
ylabel(texlabel('Power (mV^2/Hz)'))
