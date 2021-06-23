function [n_Index  fixrate] = Neuron_Data_INDEX_odrdistVar_Stim(filename)
% only calculate the 8 most used features, the other features will be
% ignored
oppClassList = [5 6 7 8 1 2 3 4 13 14 15 16 9 10 11 12];
try
    load(filename);
    try
        MatData = MUAData;
    end
    if ~isempty(MatData)
        NClass = [1:5];
        ntr = 0;
        tempRateStim =[];
        tempRate =[];
        tempRateCDStim =[];
        tempRateCD =[];
        fixRateStim =[];
        fixRate =[];
        for n = 1:length(NClass)%length(MatData.class)
            cl = NClass(n);
            if ~isempty(MatData.class(cl).ntr)
                nIndex = find([MatData.class(cl).ntr.Stim] == 1);
                tempRateStim = [tempRateStim [MatData.class(cl).ntr(nIndex).cuerate]];
                tempRateCDStim = [tempRateCDStim [MatData.class(cl).ntr(nIndex).cuedelay]];
                fixRateStim = [ fixRateStim [MatData.class(cl).ntr(nIndex).fix]];
                
                nIndex = find([MatData.class(cl).ntr.Stim] == 0);
                tempRate = [tempRate [MatData.class(cl).ntr(nIndex).cuerate]];
                tempRateCD = [tempRateCD [MatData.class(cl).ntr(nIndex).cuedelay]];
                fixRate = [fixRate [MatData.class(cl).ntr(nIndex).fix]];
                ntr = ntr +1;
            end
        end
        var1_stim = nanmean(tempRateStim);
        var1 = nanmean(tempRate);
        var1_stimcd = nanmean(tempRateCDStim);
        var1_cd = nanmean(tempRateCD);
        
        NClass = [6:10];
        ntr = 0;
        tempRateStim =[];
        tempRate =[];
        tempRateCDStim =[];
        tempRateCD =[];
%         fixRateStim =[];
%         fixRate =[];
        for n = 1:length(NClass)%length(MatData.class)
            cl = NClass(n);
            if ~isempty(MatData.class(cl).ntr)
                nIndex = find([MatData.class(cl).ntr.Stim] == 1);
                tempRateStim = [tempRateStim [MatData.class(cl).ntr(nIndex).cuerate]];
                tempRateCDStim = [tempRateCDStim [MatData.class(cl).ntr(nIndex).cuedelay]];
                fixRateStim = [ fixRateStim [MatData.class(cl).ntr(nIndex).fix]];
                
                nIndex = find([MatData.class(cl).ntr.Stim] == 0);
                tempRate = [tempRate [MatData.class(cl).ntr(nIndex).cuerate]];
                tempRateCD = [tempRateCD [MatData.class(cl).ntr(nIndex).cuedelay]];
                fixRate = [fixRate [MatData.class(cl).ntr(nIndex).fix]];
                ntr = ntr +1;
            end
        end
        var2_stim = nanmean(tempRateStim);
        var2 = nanmean(tempRate);
        var2_stimcd = nanmean(tempRateCDStim);
        var2_cd = nanmean(tempRateCD);
        varfix1_stim = nanmean(fixRateStim);
        varfix1 = nanmean(fixRate);
        %%%white task index
        n_Index(1) = abs(var1_stim-var2_stim)/(var1_stim+var2_stim);
        n_Index(2) = abs(var1_stimcd-var2_stimcd)/(var1_stimcd+var2_stimcd);
        n_Index(5) = abs(var1-var2)/(var1+var2);
        n_Index(6) = abs(var1_cd-var2_cd)/(var1_cd+var2_cd);
        
        fixrate(1) = nanmean(fixRateStim);
        fixrate(3) = nanmean(fixRate);
        %blue task
       NClass = [11:14];
        ntr = 0;
        tempRateStim =[];
        tempRate =[];
        tempRateCDStim =[];
        tempRateCD =[];
        fixRateStim =[];
        fixRate =[];
        for n = 1:length(NClass)%length(MatData.class)
            cl = NClass(n);
            if ~isempty(MatData.class(cl).ntr)
                nIndex = find([MatData.class(cl).ntr.Stim] == 1);
                tempRateStim = [tempRateStim [MatData.class(cl).ntr(nIndex).cuerate]];
                tempRateCDStim = [tempRateCDStim [MatData.class(cl).ntr(nIndex).cuedelay]];
                fixRateStim = [ fixRateStim [MatData.class(cl).ntr(nIndex).fix]];
                
                nIndex = find([MatData.class(cl).ntr.Stim] == 0);
                tempRate = [tempRate [MatData.class(cl).ntr(nIndex).cuerate]];
                tempRateCD = [tempRateCD [MatData.class(cl).ntr(nIndex).cuedelay]];
                fixRate = [fixRate [MatData.class(cl).ntr(nIndex).fix]];
                ntr = ntr +1;
            end
        end
        var1_stim = nanmean(tempRateStim);
        var1 = nanmean(tempRate);
        var1_stimcd = nanmean(tempRateCDStim);
        var1_cd = nanmean(tempRateCD);
        
        NClass = [16:19];
        ntr = 0;
        tempRateStim =[];
        tempRate =[];
        tempRateCDStim =[];
        tempRateCD =[];
%         fixRateStim =[];
%         fixRate =[];
        for n = 1:length(NClass)%length(MatData.class)
            cl = NClass(n);
            if ~isempty(MatData.class(cl).ntr)
                nIndex = find([MatData.class(cl).ntr.Stim] == 1);
                tempRateStim = [tempRateStim [MatData.class(cl).ntr(nIndex).cuerate]];
                tempRateCDStim = [tempRateCDStim [MatData.class(cl).ntr(nIndex).cuedelay]];
                fixRateStim = [ fixRateStim [MatData.class(cl).ntr(nIndex).fix]];
                
                nIndex = find([MatData.class(cl).ntr.Stim] == 0);
                tempRate = [tempRate [MatData.class(cl).ntr(nIndex).cuerate]];
                tempRateCD = [tempRateCD [MatData.class(cl).ntr(nIndex).cuedelay]];
                fixRate = [fixRate [MatData.class(cl).ntr(nIndex).fix]];
                ntr = ntr +1;
            end
        end
        var2_stim = nanmean(tempRateStim);
        var2 = nanmean(tempRate);
        var2_stimcd = nanmean(tempRateCDStim);
        var2_cd = nanmean(tempRateCD);
        varfix1_stim = nanmean(fixRateStim);
        varfix1 = nanmean(fixRate);
        
        n_Index(3) = abs(var1_stim-var2_stim)/(var1_stim+var2_stim);
        n_Index(4) = abs(var1_stimcd-var2_stimcd)/(var1_stimcd+var2_stimcd);
        n_Index(7) = abs(var1-var2)/(var1+var2);
        n_Index(8) = abs(var1_cd-var2_cd)/(var1_cd+var2_cd);
        fixrate(2) = nanmean(fixRateStim);
        fixrate(4) = nanmean(fixRate);
        
    else
        disp(['empty matdata : ' filename])
        n_Index = nan;
        n_Index2 = nan;
        max_class = nan;
        min_class = nan;fixrate = nan;
    end
catch
    disp([filename ': '  lasterr])
    n_Index = [nan nan nan nan nan nan nan nan];
    fixrate = [nan nan nan nan];
end