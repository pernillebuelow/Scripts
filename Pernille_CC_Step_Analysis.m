%% get file list
close all

%d = 'C:\Users\Michael\Desktop\New Folder'; %change this to the folder containing all the files
d = '/Volumes/PENNY/FXS data/SpikeInterval_MatLab/WT_TTXAPV/WT_TTXAPV_MatLab';

flist = ls(d);
%flist(1:2,:)=[];



flist2='';
for ind = 1:size(flist,1)
    fname = flist(ind,:);
%     while fname(end)==' '
%         fname(end)=[];
%     end
    flist2(ind,:) = ['''' fname ''' []'];
end
flist2

% flist2=[repmat('''',size(flist,1),1) flist]



%% set threshold parameters
close all

file_matrix = {
%    '032217 WT TTXAPV cl2 cell1 rec1 Export.mat' [0 0 0 0 0 0 -0.02 -0.02 -0.02 -0.02 -0.02 0.04 -0.02 -0.02 -0.02]
%     '032217 WT TTXAPV cl2 cell1 rec2 Export.mat' [-0.02]
%    % '032217 WT TTXAPV cl2 cell1 rec3 Export.mat' []
%     '032217 WT TTXAPV cl2 cell1 rec4 Export.mat' [-0.02]
%      '032217 WT TTXAPV cl2 cell2 rec1 Export.mat' [-0.02]
%      '032217 WT TTXAPV cl2 cell2 rec2 Export.mat' []
%      '032217 WT TTXAPV cl2 cell2 rec3 Export.mat' []
%      '032217 WT TTXAPV cl2 cell2 rec4 Export.mat' []
%   '040617 WT IE cL3 TTXAPV cell1 rec1 Export.mat' [-0.0214]
%     '040617 WT IE cL3 TTXAPV cell1 rec2 Export.mat' [-0.022]
%     '040617 WT IE cL3 TTXAPV cell1 rec3 Export.mat' [-0.022]
%    '040617 WT IE cL3 TTXAPV cell2 rec1 Export.mat' [0 0 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.002 -0.005 0 0 0 0 0]
%    '040617 WT IE cL3 TTXAPV cell2 rec2 Export.mat' [0 0 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.005 0 0 0 0 0]
%     '040617 WT IE cL3 TTXAPV cell2 rec3 Export.mat' [0 0 -0.01 -0.01 0 0 0 0 0 -0.001 0 0 0 0 0]
%     '040617 WT IE cL3 TTXAPV cell3 rec1 Export.mat' [-0.01]
%     '040617 WT IE cL3 TTXAPV cell3 rec2 Export.mat' [-0.01]
%     '040617 WT IE cL3 TTXAPV cell3 rec3 Export.mat' [0]
%     '051917 L8 WT cl1 TTXAPV rheo cell1 rec1 Export.mat' [0 0 -0.01 -0.01 -0.01 -0.01 -0.005 -0.005 -0.005 0 0 0 0.004 0.004 0.004]
%     '051917 L8 WT cl1 TTXAPV rheo cell1 rec2 Export.mat' [0 0 -0.01 -0.01 -0.01 -0.01 -0.005 0 0 0 0.005 0.005 0.008 0.008 0.008]
%     '051917 L8 WT cl1 TTXAPV rheo cell1 rec3 Export.mat' [0 0 -0.01 -0.01 -0.01 -0.01 0 0 0 0.005 -0.005 -0.005 0.008 0.008 0.02]
%     '051917 L8 WT cl1 TTXAPV rheo cell1 rec4 Export.mat' [0 0 -0.01 -0.01 -0.01 -0.01 0 0 0 0 0.005 -0.005 0.02 0.02 0.02]
%     '051917 L8 WT cl1 TTXAPV rheo cell2 rec1 Export.mat' [0 0 -0.01 -0.01 -0.01 -0.01 -0.005 0 -0.005 -0.004 0 0 0 0.005 0.005]
%     '051917 L8 WT cl1 TTXAPV rheo cell2 rec2 Export.mat' [0 0 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 0 0 0 0 0 0.005 0.005]
%     '051917 L8 WT cl1 TTXAPV rheo cell2 rec3 Export.mat' [0 0 -0.01 -0.01 -0.01 -0.01 -0.005 0 0 0 0 0 0.005 0.005 0.005]
%     '051917 L8 WT cl1 TTXAPV rheo cell3 rec1 Export.mat' [0 0 -0.01 -0.01 -0.005 0 -0.004 0 0.02 0.02 0.02 0.02 0.02 0.02 0.02]
%     '051917 L8 WT cl1 TTXAPV rheo cell3 rec2 Export.mat' [0 0 -0.01 -0.01 -0.01 -0.005 -0.005 0 0.02 0.02 0.02 0.02 0.02 0.02 0.02]
%     '051917 L8 WT cl1 TTXAPV rheo cell3 rec3 Export.mat' [0 0 -0.01 -0.01 -0.01 -0.005 -0.005 0 0.02 0.02 0.02 0.02 0.02 0.02 0.02]
%     '051917 L8 WT cl4 TTXAPV rheo cell2 rec1 Export.mat' [0 0 -0.01 -0.01 -0.01 -0.01 -0.01 -0.005 -0.005 -0.005 -0.005 -0.005 -0.005 -0.005 -0.003]
%     '051917 L8 WT cl4 TTXAPV rheo cell2 rec2 Export.mat' [0 0 -0.01 -0.01 -0.01 -0.01 -0.01 -0.005 -0.005 -0.005 -0.005 -0.005 -0.005 -0.005 -0.005]
%     '051917 L8 WT cl4 TTXAPV rheo cell2 rec3 Export.mat' [0 0 -0.015 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.005 -0.005 -0.005 -0.005 -0.005 -0.005]
%    '051917 L8 WT cl4 TTXAPV rheo cell2 rec4 Export.mat' [0 0 -0.01 -0.015 -0.01 -0.01 -0.01 -0.01 -0.005 -0.005 -0.005 -0.005 -0.005 -0.005 -0.005]
%     '051917 L8 WT cl4 TTXAPV rheo cell3 rec1 Export.mat' [0 0.03 0.03 0.03 -0.01 -0.005 -0.005 0 0 0.005 0.005 0.005 0.005 0.005 0.005]
%     '051917 L8 WT cl4 TTXAPV rheo cell3 rec2 Export.mat' [0 -0.01 -0.01 -0.01 0 0 0 0 0 0 0 0.005 0.005 0.005 0.01]
%     '051917 L8 WT cl4 TTXAPV rheo cell3 rec3 Export.mat' [0 -0.01 -0.015 0 0 0 0 0 0 0.005 0.005 0.02 0.02 0.02 0.02]
%     '070717_L9 WT cl4-TTXAPV_c1-rec1  Export.mat' [0 -0.01 -0.015 -0.015 -0.015 -0.015 -0.015 -0.015 -0.015 -0.0129 -0.015 -0.015 -0.013 -0.013 -0.014]
%     '070717_L9 WT cl4-TTXAPV_c1-rec2  Export.mat' [-0.0146]  
%     '070717_L9 WT cl4-TTXAPV_c1-rec3  Export.mat' [-0.015]
%     '070717_L9 WT cl4-TTXAPV_c1-rec4  Export.mat' [-0.015]
%     '070717_L9 WT cl4-TTXAPV_c2-rec1  Export.mat' [0 0 0 0 -0.02 -0.02 -0.02 -0.02 -0.017 -0.02 -0.015 -0.015 0 0 0]
%     '070717_L9 WT cl4-TTXAPV_c2-rec2  Export.mat' [0 0 0 0 -0.02 -0.025 -0.02 -0.02 -0.018 -0.015 0 0 0 0 0]
%     '070717_L9 WT cl4-TTXAPV_c2-rec3 Export.mat' [0 0 0 0 -0.025 -0.02 -0.02 -0.02 -0.02 -0.02 0 0 0 0 0]
%      '070717_L9 WT cl4-TTXAPV_c3-rec1  Export.mat' [0 0 0 -0.01 -0.01 -0.015 -0.005 -0.005 0 0 0 0 0 0 0]
%     '070717_L9 WT cl4-TTXAPV_c3-rec2  Export.mat' [0 0 0 -0.015 -0.015 -0.015 -0.015 -0.01 -0.01 -0.01 -0.008 -0.005 -0.005 -0.005 -0.005]
%     '070717_L9 WT cl4-TTXAPV_c3-rec3  Export.mat' [0 0 0 0 0 -0.015 -0.015 -0.015 -0.015 -0.01 0 0 0 0 0]
%     '071217_IE_WT_cl2-TTXAPV-cell1-rec1 Export.mat' [0 0 0 0 0 0 0 0 0 -0.01 -0.01 -0.01 -0.01 -0.01 0]
%     '071217_IE_WT_cl2-TTXAPV-cell1-rec2 Export.mat' [0 0 0 0 0 -0.01 -0.01 -0.02 -0.01 -0.01 -0.01 -0.005 -0.005 -0.005 -0.005]
%     '071217_IE_WT_cl2-TTXAPV-cell1-rec3 Export.mat' [0 0 0 0 0 -0.01 -0.01 -0.02 -0.01 -0.01 -0.01 -0.005 -0.005 -0.005 -0.005]
%     '071217_IE_WT_cl2-TTXAPV-cell1-rec4 Export.mat' [0 0 0 0 0 -0.02 -0.015 -0.015 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.005]
%      '071217_IE_WT_cl2-TTXAPV-cell2-rec1 Export.mat' [-0.01]
%     '071217_IE_WT_cl2-TTXAPV-cell2-rec2 Export.mat' [-0.01]
%     '071217_IE_WT_cl2-TTXAPV-cell2-rec3 Export.mat' [-0.01]
%     '071217_IE_WT_cl2-TTXAPV-cell3-rec1 Export.mat' [-0.015]
%      '071217_IE_WT_cl2-TTXAPV-cell3-rec2 Export.mat' []
%     '071217_IE_WT_cl2-TTXAPV-cell3-rec3 Export.mat' []
%     '101917_L16_IE_WT_cl2-TTXAPV-cell1-rec1 Export.mat' [0 0 0 -0.01 -0.013 -0.01 -0.011 -0.01 -0.01 -0.01 0 0 0 0 0]
%     '101917_L16_IE_WT_cl2-TTXAPV-cell1-rec3 Export.mat' [0 0 0 -0.01 -0.02 -0.01 -0.01 -0.015 -0.015 -0.015 0 0 0 0 0]
%      '101917_L16_IE_WT_cl2-TTXAPV-cell1-rec5 Export.mat' [0 0 0 0 -0.01 -0.02 -0.02 -0.02 0 0 0 0 0 0 0]
%      '101917_L16_IE_WT_cl2-TTXAPV-cell2-rec1 Export.mat' [0 0 0 0 -0.01 -0.015 -0.02 -0.02 -0.02 -0.015 -0.015 -0.017 -0.017 -0.017 -0.017]
%      '101917_L16_IE_WT_cl2-TTXAPV-cell2-rec2 Export.mat' [0 0 0 0 -0.01 -0.015 -0.02 -0.02 -0.022 -0.02 -0.017 -0.016 -0.016 -0.016 -0.016]
%      '101917_L16_IE_WT_cl2-TTXAPV-cell2-rec3 Export.mat' [0 0 0 0 -0.01 -0.015 -0.015 -0.015 -0.015 -0.015 -0.015 -0.015 -0.015 -0.015 -0.015]
%      '101917_L16_IE_WT_cl2-TTXAPV-cell3-rec1 Export.mat' [-0.02]
%      '101917_L16_IE_WT_cl2-TTXAPV-cell3-rec5 Export.mat' [-0.02]
%      '101917_L16_IE_WT_cl2-TTXAPV-cell3-rec6 Export.mat' [-0.02]
%      '101917_L16_IE_WT_cl2-TTXAPV-cell4-rec1 Export.mat' [-0.006]
%      '101917_L16_IE_WT_cl2-TTXAPV-cell4-rec2 Export.mat' [-0.01]
%      '101917_L16_IE_WT_cl2-TTXAPV-cell4-rec3 Export.mat' [-0.01]
%      '101917_L16_IE_WT_cl2-TTXAPV-cell4-rec4 Export.mat' [-0.01]
%      '101917_L16_IE_WT_cl2-TTXAPV-cell5-rec1 Export.mat' [-0.02]
%      '101917_L16_IE_WT_cl2-TTXAPV-cell5-rec2 Export.mat' [-0.02]
%      '101917_L16_IE_WT_cl2-TTXAPV-cell5-rec3 Export.mat' [-0.02]
%      '4465b Export.mat' [0 0.04 0 0 0 -0.01 -0.01 -0.01 0.04 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01]
%      '4465c Export.mat' [-0.01]
     '725b Export.mat' [0.04 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
     '725c Export.mat' [0 0 0 0 0 0.04 0 0 0 0 0.04 0 0 0 0]
     '725d Export.mat' [0.04 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
    




};


    


type = 'n';

% for each file
FR_mean_tot = [];
FR_max_tot = [];
peak_V = nan(50,15,size(file_matrix,1));
AHP_V = peak_V;

for fnd = 1:size(file_matrix,1) 
    fname = file_matrix{fnd,1};
    while fname(end)==' '
        fname(end)=[];
    end
    load([d '/' fname])

    thresh_params = file_matrix{fnd,2};
    
% load 'C:\Users\Michael\Desktop\894b Export.mat'
% load 'C:\Users\Michael\Desktop\932b Exportx.mat'


% read in file
t = c001_Time;
V = [];
I = [];
si = (t(2)-t(1));

for ind = 1:15
    V_ID = num2str(ind*2);
    if length(V_ID)==1,V_ID = ['0' V_ID];end
    I_ID = num2str(ind*2+1);
    if length(I_ID)==1,I_ID = ['0' I_ID];end
    
    eval(['V(:,ind) = c0' V_ID '_Voltage;'])
    eval(['I(:,ind) = c0' I_ID '_Current;'])
end

FR_max = zeros(1,size(V,2));
FR_mean = FR_max;
adapt_ratio = FR_max;
slope_max = FR_max;
AUC = FR_max;
tail = FR_max;
burst_dur = FR_max;

%for each current step
for jnd = 1:size(V,2)
    trace = V(:,jnd);
   if t(end)==1
            base_index = find(t==0.3);
            start_index = find(t==0.4);
            end_index = find(t==0.9);
            baseline = mean(trace(base_index:start_index));
            auc_index = start_index:end_index;
            tail_index = end_index:length(trace);
   else
            start_index = find(t==0.1);
            end_index = find(t==0.6);
            baseline = mean(trace(1:start_index));
            auc_index = start_index:end_index;
            tail_index = end_index:length(trace);
   end
   AUC(jnd) = sum(trace(auc_index)-baseline)*si;
   tail(jnd) = sum(trace(tail_index)-baseline)*si;
   
    if isempty(thresh_params)
        thresh = 0.8*max(V(:)) + 0.2*min(V(:));   
    elseif length(thresh_params)==1
        thresh = thresh_params;
    elseif length(thresh_params)==15
        thresh = thresh_params(jnd);
    end
%     thresh=-.05;

    ds=find(diff(trace>thresh)==1);
    de=find(diff(trace>thresh)==-1);
    peaks=[];
    AHP = [];
    AHP_index = [];
    for ind = 1:length(de)
        [mx mi] = max(trace(ds(ind):de(ind)));
    %     if ds(ind)>=stepstart && de(ind)<=stepend && mx>-10
            peaks(ind)=ds(ind)+mi-1;
    %     end
    end
    peaks(peaks==0)=[];
    n_spikes(jnd) = length(peaks);
    IFR = 1./diff(t(peaks));
    
%     if jnd==10   %change this to look at different steps
%         figure(1000), hold on
%         plot(t(peaks(1:end-1)),IFR)
%     end
        
    peak_V(1:length(peaks),jnd,fnd)=trace(peaks);
    
    for ind = 1:(length(peaks)-1)
        [mn mi]= min(trace(peaks(ind):peaks(ind+1)));
        AHP(ind)=mn;
        AHP_index(ind) = mi+peaks(ind)-1;
    end
    
    peak_V(1:length(peaks),jnd,fnd)=trace(peaks);
    AHP_V(1:length(AHP),jnd,fnd)=AHP;

    if ~isempty (IFR)
        FR_max(jnd) = IFR(1);
        FR_mean(jnd) = mean(IFR);
        adapt_ratio(jnd) = IFR(1)/IFR(end);
        burst_dur(jnd) = (peaks(end)-peaks(1))*si;
    end
    if ~isempty(peaks)
        [sm smi] = max(diff(trace(start_index:peaks(1))));
        slope_max(jnd) = sm/si; %V/s (mV/ms)
        plot((smi+start_index)*si+(jnd-1),trace(smi+start_index),'mo')
    end
        
    figure(fnd)
    plot(t+(jnd-1),trace),hold on, plot(t(peaks)+(jnd-1),trace(peaks),'ro')
    fill([t(auc_index) t(auc_index(end))]+(jnd-1),[trace(auc_index); baseline],'b')
    fill([t(tail_index(1)) t(tail_index) t(tail_index(end))]+(jnd-1), [baseline; trace(tail_index); baseline],'r')  
    plot(t(AHP_index)+(jnd-1),trace(AHP_index),'mx')
    
    %plot(t(2:end),diff(trace)),hold on
    title(fname)
 
    %figure(2),plot(t(peaks(1:(end-1))),IFR),hold on

end

FR_max(FR_max==0)=nan;
FR_mean(FR_mean==0)=nan;

% figure
% plot(V(:,9),'r'),hold on
% plot(V(:,15))
% figure, plot(FR_max), hold on, plot(FR_mean), plot(IFR)

FR_mean_tot(:,fnd)=FR_mean';
FR_max_tot(:,fnd)=FR_max';
adapt_ratio_tot(:,fnd)=adapt_ratio';
slope_max_tot(:,fnd)=slope_max';
AUC_tot(:,fnd)=AUC;
tail_tot(:,fnd)=tail;
n_spikes_tot(:,fnd)=n_spikes;
burst_dur_tot(:,fnd)=burst_dur;
end


% for ind = 1:size(FR_mean_tot,1)
%     vec = FR_mean_tot(ind,:);
%     vec(isnan(vec))=[];
%     FR_mean_tot_mean(ind)=mean(vec)
% end

FR_mean_tot
FR_max_tot
peak_V = reshape(peak_V,50,[]);
AHP_V = reshape(AHP_V,50,[]);
