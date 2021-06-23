function [psthM psthNM fixMean psthCT1 psthCT2 isReverse maxPsth maxCpsth nnm] = extractPSTH_OdrdistVar_StimOneFile(filename,bin_edges,bin_width,IsSmooth,phase,varargin)
% divide class by number of 1st stimuli
%smooth parameters
% Sim and non-Stim trials in one file
%seperate them to two output
start_t= -1;
end_t= 6;
smooth_sdv = 0.080;

load(filename)
% load(filename(1:8));% behavior data
fixationTime  = 1;
cuerate = [];
fixrate =[];
allTSm = [];
allTSnm = [];
nnm = 0;
nnn = 0;
% ClassStructure = AllData.ClassStructure;
nn =0;
Nclass = length(MatData.class);
maxPsth =nan;
% reverseRF = 0;
for n = 1: Nclass
    allTSm(n).TS =[];
    allTSnm(n).TS =[];
    nnm(n) =0;
    nnn(n) =0;
end
%%%% find the class psth will be extracted from
for m = 1 :10;% length(MatData.class)
    %     cuerate(m) = nanmean([MatData.class(m).ntr.cuedelay]);
    if ~isempty(MatData.class(m).ntr)
    fixrate = [fixrate [MatData.class(m).ntr.fix]];
    end
end
fixMean(1,1) = nanmean(fixrate);
fixrate =[];

for m = 11 :20;% length(MatData.class)
    %     cuerate(m) = nanmean([MatData.class(m).ntr.cuedelay]);
    if ~isempty(MatData.class(m).ntr)
    fixrate = [fixrate [MatData.class(m).ntr.fix]];
    end
end
fixMean(1,2) = nanmean(fixrate);

cuedelayT1=[];
for m = [1:4]
     if ~isempty(MatData.class(m).ntr)
    eval(['cuedelayT1 =[cuedelayT1 [MatData.class(m).ntr.' phase ']];'])
     end
end
cuedelayT2=[];
for m = [6:9]
     if ~isempty(MatData.class(m).ntr)
    eval(['cuedelayT2 =[cuedelayT2 [MatData.class(m).ntr.'  phase ']];'])
     end
end

cuedelayT3=[];
for m = [11:14]
     if ~isempty(MatData.class(m).ntr)
    eval(['cuedelayT3 =[cuedelayT3 [MatData.class(m).ntr.' phase ']];'])
     end
end
cuedelayT4=[];
for m = [16:19]
     if ~isempty(MatData.class(m).ntr)
    eval(['cuedelayT4 =[cuedelayT4 [MatData.class(m).ntr.' phase ']];'])
     end
end

[a b] = max([ nanmean(cuedelayT1)  nanmean(cuedelayT2)  nanmean(cuedelayT3)  nanmean(cuedelayT4)]);
% if nanmean(cuedelayT1) < nanmean(cuedelayT2)
if mod(b,2)==0;
    isReverse ='r';
else
    isReverse ='n';
end
% [maxrate maxClass] = max(cuerate(1:2));%decide RF
% if maxClass >1
%     reverseRF = 'r';
% else
%     reverseRF = 'n;;
% end
if strcmpi(isReverse,'r')
    % if isReverse
    reverseN = [6 7 8 9 10 1 2 3 4 5 16 17 18 19 20 11 12 13 14 15];
end
%%%% Main
try
    for m = 1 : length(MatData.class)
        if ~isempty(MatData.class(m).ntr)
            for nn = 1:length(MatData.class(m).ntr)
                if ~MatData.class(m).ntr(nn).Stim
                    TS = MatData.class(m).ntr(nn).TS-(MatData.class(m).ntr(nn).Cue_onT-fixationTime);
                    allTSm(m).TS = [allTSm(m).TS  TS];
                    nnm(m) = nnm(m)+1;
                else
                    TS = MatData.class(m).ntr(nn).TS-(MatData.class(m).ntr(nn).Cue_onT-fixationTime);
                    allTSnm(m).TS = [allTSnm(m).TS  TS];
                    nnn(m) = nnn(m)+1;
                end
            end
        else
            allTSm(m).TS = [];
             allTSnm(m).TS =[];
        end
    end
catch
    lasterr
    psthM(m).psth =  ones(1,size(bin_edges,2)) *nan;
    psthNM(m).psth =  ones(1,size(bin_edges,2)) *nan;
    psthCT1=ones(1,size(bin_edges,2)) *nan;
    psthCT2=ones(1,size(bin_edges,2)) *nan;
    maxPsth =nan;

    %     varargout{1} =psthNM;
    return
end
maxPsth=0;
for m = 1 : Nclass
    if nnn(m) >= 1%~isempty(allTSm(m).TS)
        if IsSmooth
            psthM(m).psth = sdfpsth(allTSm(m).TS,start_t,end_t+bin_width,bin_width,smooth_sdv)/(nnm(m)*bin_width);
            psthNM(m).psth = sdfpsth(allTSnm(m).TS,start_t,end_t+bin_width,bin_width,smooth_sdv)/(nnn(m)*bin_width);
        else
            psthM(m).psth  = histc(allTSm(m).TS,bin_edges)/(bin_width*(nnm(m)));
            psthNM(m).psth  = histc(allTSnm(m).TS,bin_edges)/(bin_width*(nnn(m)));
        end
        maxPsth = max([max(psthNM(m).psth) ,maxPsth,max(psthM(m).psth) ]);
    else
        psthM(m).psth  = ones(1,size(bin_edges,2)) * nan;
        psthNM(m).psth  = ones(1,size(bin_edges,2)) * nan;
%         psthCT1=ones(size(bin_edges,2),1) *nan;
%         psthCT2=ones(size(bin_edges,2),1) *nan;

    end
    
end
if  strcmpi(isReverse,'r')
    psthM = psthM(reverseN);
    psthNM = psthNM(reverseN);
    nnm=nnm(reverseN);
end
% varargout{1} =psthNM;

classList = {[1:4]; [6:9]; [11:14]; [16:19]};
maxCpsth =0;
for n = 1 : length(classList)
    allTSC1 =[];ntr1 =0;%stim
    allTSC2 =[];ntr2 =0;%no-stim
    classes = classList{n};
    psth1 =[];
     psth2 =[];
    for m = classes
        psth1 =[psth1; psthM(m).psth];
         psth2 =[psth2; psthNM(m).psth];
%         allTSC1 = [allTSC1 allTSm(m).TS];
%         ntr1= ntr1 +  nnm(m);
%         allTSC2 = [allTSC2 allTSnm(m).TS];
%         ntr2= ntr2 +  nnn(m);
    end
psthCT1(n).psth =nanmean( psth1);
psthCT2(n).psth =nanmean( psth2);
end
% for n = 1 : length(classList)
%     allTSC1 =[];ntr1 =0;%stim
%     allTSC2 =[];ntr2 =0;%no-stim
%     classes = classList{n};
%     for m = classes
%         allTSC1 = [allTSC1 allTSm(m).TS];
%         ntr1= ntr1 +  nnm(m);
%         allTSC2 = [allTSC2 allTSnm(m).TS];
%         ntr2= ntr2 +  nnn(m);
%     end
%     if ~isempty(allTSC1)
%         if IsSmooth
%             psthCT1(n).psth = sdfpsth(allTSC1,start_t,end_t+bin_width,bin_width,smooth_sdv)/(ntr1*bin_width);
%             psthCT2(n).psth = sdfpsth(allTSC2,start_t,end_t+bin_width,bin_width,smooth_sdv)/(ntr2*bin_width);
%         else
%             psthCT1(n).psth  = histc(allTSC1,bin_edges)/(bin_width*ntr1);
%             psthCT2(n).psth  = histc(allTSC2,bin_edges)/(bin_width*ntr2);
%         end
%         maxCpsth=max([maxCpsth max(psthCT2(n).psth)  max(psthCT1(n).psth)]);
%     else
%         psthCT1(n).psth  = ones(size(bin_edges,2),1) * nan;
%         psthCT2(n).psth  = ones(size(bin_edges,2),1) * nan;
%     end
% end
% 
% if  strcmpi(isReverse,'r')
%     psthCT1 = psthCT1([2 1 4 3]);
%     psthCT2 = psthCT2([2 1 4 3]);
% end