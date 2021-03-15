function [SPK_trg,PSTH,var_idx,chirp_all, bar_all, color_all,flash_all,data,nr_name,Bad_cell_name,e_Sta_recP]=averaging_func(ds_list)


% claculate the average of trials and make a matrix for each stimulus. The input is the location of
% matlab spike times file and the output is a matrix of normalized psth of
% all data.


% outputs:

%SPK_trg:       spike times of chirp data
%PSTH:          PSTH of chirp data binsize 100 ms
%var_idx:       variance ratio(explained in paper)
%chirp_all:     median of normalized pseudoclcum traces for chirp data (normalised to max(abs()))
%bar_all:       average of normalized pseudoclcum traces for movingbar data (normalised to max(abs()))
%color_all      average of normalized pseudoclcum traces for color data (normalised to max(abs()))
%flash_all      average of normalized pseudoclcum traces for flash data (normalised to max(abs()))
%data           moving bar data average and by repetition with nearby electrodes data
%nr_name        recording data and channel name of each unit
%Bad_cell_name: cells with no activity 
%e_Sta_recP:    electrical STAs of each unit

%%


%pst='psth';% if you want to have nomal psth
pst='psth_ogb';% if you want to plot ogb psth
m_bar=[];
chirp_all=[];flash_all=[];
unnormalized_chirp=[];unnormalized_flash=[];unnormalized_bar=[];unnormalized_color=[];
bar_all=[];color_all=[];fig_nmbr=0;cnt=0;

nrn_cnt=0;
trg_data=[];trg_label2=[];trg_label=[];
c=0;
i22=0;
var_idx=[];
e_Sta_rec={};

for ds=1:size(ds_list,2)
    %% compute m-bar  average
    chcnt=0;
    bar=load([ds_list(ds).FolderName,'Matlab_results','\Moving Bar.mat']);% load m-bar
    
    staname=['\e_STA ' ,ds_list(ds).FileName];
    sta=load([ds_list(ds).FolderName,'Matlab_results',staname]);%load STA
    
    mydir = [ds_list(ds).FolderName,'Matlab_results','\Moving Bar.mat'];% load HT
    idcs = strfind (mydir, filesep);
    newdir = mydir (1: idcs(end-1)-1);
    foldername = newdir(idcs(end-2)+1:end) ;
    load([newdir ,'\HT_',foldername]);
    foldername = regexprep(foldername, ' ', '_');
    
brn_cnt=0;
    for i2=1:size(bar.Data,2)-1
        
        
        if length(bar.Data(i2).spks)>1 
        i2;cnt=cnt+1;
        e_Sta_recP(cnt).esta=sta.Data(i2).e_sta';
        e_Sta_recP(cnt).name=sta.Data(i2).name';
        e_Sta_recP(cnt).FolderName=ds_list(ds).FolderName;

%        aa=[ds cnt];
        ch=char(bar.Data(i2).name());
        ch_name=ch(1:end-1);
        nrn_cnt=nrn_cnt+1; brn_cnt=brn_cnt+1;

        nr_name{nrn_cnt,1}=[foldername, ch];
        bdnr_name{brn_cnt,1}=[foldername, ch];

        for i=1:8% orientations
            deg=fieldnames(bar.Data(i2).psth);
            deg_name=deg{i};
            nrst_stm=HT.(ch_name).(deg_name).nrst_stm;

            nrst_stm_and_neighbors=[nrst_stm-1 nrst_stm nrst_stm+1];% add the data of neighboring stimuli
            if nrst_stm==1
                nrst_stm_and_neighbors(1)=[];
            end

            if nrst_stm>=length(find(~cellfun(@isempty,{bar.Data(i2).(pst).(deg_name)})))
                nrst_stm_and_neighbors(end)=[];
            end


            rmov_bar_8_ogb= bar.Data(i2).(pst)(nrst_stm).(deg_name);
            rmov_bar_8_ogb(isnan(rmov_bar_8_ogb))=0;
            mov_bar_8_ogb(i,:)=rmov_bar_8_ogb;%/max(abs(rmov_bar_8_ogb));%normalize
            mov_bar_1_plus_neighbor_ogb_trls=[];mov_bar_1_plus_neighbor_ogb=[];
            mov_bar_1_plus_neighbor=[];
            for nb=1:length(nrst_stm_and_neighbors)% considering neighboring stimuli
                rd= bar.Data(i2).(pst)(nrst_stm_and_neighbors(nb)).(deg_name);
                rd(isnan(rd))=0;
                
                rd_trl= bar.Data(i2).([pst,'_trls'])(nrst_stm_and_neighbors(nb)).(deg_name);

                %nd=rd/max(abs(rd));%normalize
                mov_bar_1_plus_neighbor_ogb(nb,:)= rd;
                if sum(sum(rd_trl(sum(rd_trl,2)>0,:)))>1
                mov_bar_1_plus_neighbor_ogb_trls(:,:,nb)= rd_trl;%(sum(rd_trl,2)>0,:);
                end
                mov_bar_1_plus_neighbor=[rd_trl;mov_bar_1_plus_neighbor];
            end
            mov_bar_1_plus_neighbor_8(i,:)=median(mov_bar_1_plus_neighbor,1);            
            if isempty(mov_bar_1_plus_neighbor_ogb_trls)
                mov_bar_1_plus_neighbor_ogb_trls=zeros(size(rd_trl));
            end
            mov_bar_8_plus_neighbor_ogb_trls{i}=nanmedian(mov_bar_1_plus_neighbor_ogb_trls,3);
            
            mov_bar_8_plus_neighbor_ogb(i,:)= nanmean(mov_bar_1_plus_neighbor_ogb(nansum(mov_bar_1_plus_neighbor_ogb,2)>0,:),1);% average over nonzero trials. Put Nan if all trials are zero
            mov_bar_8_plus_neighbor_ogb(isnan(mov_bar_8_plus_neighbor_ogb))=0;
            mbar_nbrs=nanmean(mov_bar_1_plus_neighbor_ogb,1);
            nmbar_nbrs=mbar_nbrs/max(abs(mbar_nbrs));
            nmbar_nbrs(isnan(nmbar_nbrs))=0;
            mov_bar_8_plus_neighbor_ogb_nonan(i,:)=nmbar_nbrs;

        end
        avg_mov_bar=nanmean(mov_bar_8_ogb,1);
        data.(['data_',foldername])(i2).name=bar.Data(i2).name;
        data.(['data_',foldername])(i2).mov_bar_avg_ogb= nanmean(mov_bar_8_ogb(sum(mov_bar_8_ogb,2)>0,:),1);
        data.(['data_',foldername])(i2).mov_bar_avg_plus_neighbor_ogb= nanmean(mov_bar_8_plus_neighbor_ogb(sum(mov_bar_8_plus_neighbor_ogb,2)>0,:),1);
        data.(['data_',foldername])(i2).mov_bar_all_ogb= mov_bar_8_ogb;
        data.(['data_',foldername])(i2).mov_bar_all_plus_neighbor_ogb=mov_bar_1_plus_neighbor_8 ;%mov_bar_8_plus_neighbor_ogb;%(nansum(mov_bar_8_plus_neighbor_ogb,2)>0,:);
        data.(['data_',foldername])(i2).mov_bar_all_plus_neighbor_ogb_trls= mov_bar_8_plus_neighbor_ogb_trls;

        m_bar(i2,:)= nanmean(mov_bar_8_plus_neighbor_ogb(nansum(mov_bar_8_plus_neighbor_ogb,2)>0,:),1);% each row one channl
        m_bar_8= mov_bar_8_plus_neighbor_ogb_nonan;
        m_bar_8(end+1,:)=0;
        % draw polar plots
        if sum(m_bar_8(:))
            chcnt=chcnt+1;
           if chcnt>20% number of plots in each figure
               fig_nmbr=fig_nmbr+1;
               chcnt=1;
           end
             

%             hfig=figure(ds+fig_nmbr); %if you want to see polarplots
%             subplot(4,5,chcnt)
%             R = linspace(min(m_bar_8(:)),max(m_bar_8(:)),size(m_bar_8,2));
%             theta = linspace(0,360,9);
%             polarPcolor(double(R),theta,double(m_bar_8),'Ncircles',2,'colBar',0,'Nspokes',9,'labelR','     FR','RtickLabel',{'',''})
% 
%             title(ch_name(end-1:end))

        end
    end%% m-bar  average
    end
    %print(hfig,'-dpng',[fpath,'Bsline removed Moving Bar polar plot ',mbar_loc{ds},'.png'],'-r300')
    empty_bar=find(mean(isnan(m_bar),2));% index of cells with no activity for moving bar stim
    bad_cells=[];
    chirp=load([ds_list(ds).FolderName,'Matlab_results','\Chirp.mat']);% load chirp
    
    empty_chirp = find(cell2mat(arrayfun(@(x) isempty(x.(pst)),chirp.Data,'UniformOutput',false)));
    empty_chirp(end)=[];
    l_chirp = cell2mat(arrayfun(@(x) length(x.(pst)),chirp.Data,'UniformOutput',false));
    %neg_chirp = cell2mat(arrayfun(@(x) sum(x.(pst)),chirp.Data,'UniformOutput',false));

    %neg_chirp=neg_chirp(l_chirp~=0);
    l_chirp=l_chirp(l_chirp~=0);
    
%bd_lbl=find(neg_chirp==-l_chirp);
%% spikes of chirp
stim_dur=32;
trgrs=chirp.Data(end).spks.Chirp  ;
between_trials =find(diff(trgrs)>2*stim_dur);% detect start of different recording session in one experiment using the timing between trials
trgs=[trgrs(1:end-1),trgrs(2:end)];% modify triggers timing
trgs(end+1,:)=[trgrs(end),trgrs(end)+stim_dur];% add the last trg stim_dur=4 sec for flash

%trgs=trgs-stim_dur/2%  corrcttion for falling edge detection problem for
%flash stimulus


trgs(between_trials,2)= trgrs(between_trials,1)+stim_dur;% correct the wrong value of last trigger of each session

    baseline_time=1   ;
    chirp_avg=[];
    for i2=1:size(chirp.Data,2)-1% chirp avg
        i2;
        if isempty(chirp.Data(i2).(pst))
            chirp_avg(i2,:)=zeros(1,l_chirp(1));
        end
        chirp_avg=[chirp_avg;chirp.Data(i2).(pst)];
        
        
    ndata=chirp.Data(i2).spks;
    psth=[];psth_bs=[];bwid=.1;
    for trg=1:length(trgrs)
        trg;
        %trg_data_baselie=[trg_data_baselie ;ndata(ndata>trgs(trg,1)-baseline_time & ndata<trgs(trg,1))-trgs(trg,1)+baseline_time];
        trg_data=[trg_data ;ndata(ndata>(trgs(trg,1)) & ndata<trgs(trg,2))-trgs(trg,1)];
        trg_label=[trg_label ;i22+ones(length(ndata(ndata>trgs(trg,1) & ndata<trgs(trg,2))),1)];
        trg_label2=[trg_label2 ;c+ones(length(ndata(ndata>(trgs(trg,1)) & ndata<trgs(trg,2))),1)];
        c=c+1;% nuber of trg of all channels

            psth=[psth; (1/bwid)*histcounts(ndata(ndata>trgs(trg,1) & ndata<trgs(trg,2)),[trgs(trg,1): bwid:trgs(trg,1)+stim_dur])];
            
            psth_bs=[psth_bs; (1/bwid)*histcounts(ndata(ndata>trgs(trg,1)-baseline_time & ndata<trgs(trg,1))...
                ,[trgs(trg,1)-baseline_time: bwid:trgs(trg,1)])];

    end
       i22=1+i22;

    rmv=find(sum(psth,2)==0);% remove no data trials
    psth_bs(rmv,:)=[];
    psth(rmv,:)=[];
   %var_idx=[var_idx; nanmedian(var(psth_bs')./var(psth'))];
   %var_idx=[var_idx; median(var(psth'))/nanmedian(var([psth_bs psth]'))];
   var_idx=[var_idx;mean(var(psth))/var(psth(:))];% ratio of whole var to each time segment var
   PSTH(i22,:)=mean(psth);    
    end
    

    
    
    color=load([ds_list(ds).FolderName,'Matlab_results','\BG.mat']);% load color data
    empty_color = find(cell2mat(arrayfun(@(x) isempty(x.(pst)),color.Data,'UniformOutput',false)));
    empty_color(end)=[];

    l_color = cell2mat(arrayfun(@(x) length(x.(pst)),color.Data,'UniformOutput',false));
    l_color=l_color(l_color~=0);
    color_avg=[];
    for i2=1:size(color.Data,2)-1% color avg
        if isempty(color.Data(i2).(pst))
            color_avg(i2,:)=zeros(1,l_color(1));
        end
        color_avg=[color_avg;color.Data(i2).(pst)];
    end
    
    flash=load([ds_list(ds).FolderName,'Matlab_results','\Flash.mat']);% load flash data
    empty_flash = find(cell2mat(arrayfun(@(x) isempty(x.(pst)),flash.Data,'UniformOutput',false)));
    empty_flash(end)=[];
    l_flash = cell2mat(arrayfun(@(x) length(x.(pst)),flash.Data,'UniformOutput',false));
    l_flash=l_flash(l_flash~=0);
    flash_avg=[];
    %% spikes of flash
    stim_dur=4;
% trgrs=flash.Data(end).spks.Flash  ;
% between_trials =find(diff(trgrs)>2*stim_dur);% detect start of different recording session in one experiment using the timing between trials
% trgs=[trgrs(1:end-1),trgrs(2:end)];% modify triggers timing
% trgs(end+1,:)=[trgrs(end),trgrs(end)+stim_dur];% add the last trg stim_dur=4 sec for flash
% 
% %trgs=trgs-stim_dur/2%  corrcttion for falling edge detection problem for
% %flash stimulus
% 
% 
% trgs(between_trials,2)= trgrs(between_trials,1)+stim_dur;% correct the wrong value of last trigger of each session

  %  baseline_time=1   ;
    for i2=1:size(flash.Data,2)-1% flash avg
        if isempty(flash.Data(i2).(pst))
            flash_avg(i2,:)=zeros(1,l_flash(1));
        end
        flash_avg=[flash_avg;flash.Data(i2).(pst)];
%     ndata=flash.Data(i2).spks;
% psth=[];psth_bs=[];bwid=.1;
%     for trg=1:length(trgrs)
%         trg;
%         trg_data=[trg_data ;ndata(ndata>(trgs(trg,1)) & ndata<trgs(trg,2))-trgs(trg,1)];
%         trg_label2=[trg_label2 ;c+ones(length(ndata(ndata>(trgs(trg,1)-baseline_time) & ndata<trgs(trg,2))),1)];
%         c=c+1;% nuber of trg of all channels
%             psth=[psth; (1/bwid)*histcounts(ndata(ndata>trgs(trg,1) & ndata<trgs(trg,2)),[trgs(trg,1): bwid:trgs(trg,1)+stim_dur])];
%             
%             psth_bs=[psth_bs; (1/bwid)*histcounts(ndata(ndata>trgs(trg,1)-baseline_time & ndata<trgs(trg,1))...
%                 ,[trgs(trg,1)-baseline_time: bwid:trgs(trg,1)])];
% 
%     end
%     rmv=find(sum(psth,2)==0);% remove no data trials
%     psth_bs(rmv,:)=[];
%     psth(rmv,:)=[];
   %var_idx=[var_idx; nanmedian(var(psth_bs')./var(psth'))];
%   var_idx=[var_idx; median(var(psth'))/nanmedian(var([psth_bs psth]'))];
    end
    bad_cells=unique([empty_chirp';empty_color';empty_flash']);% all cells with no activity 
     
    bad_cells=unique([reshape(empty_bar,[], 1) ;reshape(empty_color,[], 1);reshape(empty_chirp,[], 1)]);% all cells with no activity 
    Bad_cell_name{ds}=bdnr_name(bad_cells);
    normalized_bar=[(m_bar)'./max(abs(m_bar)')]';
    normalized_chirp=[(chirp_avg)'./max(abs(chirp_avg)')]';
    normalized_color=[(color_avg)'./max(abs(color_avg)')]';
    normalized_flash=[(flash_avg)'./max(abs(flash_avg)')]';

 
    
    m_bar=[];chirp_avg=[];color_avg=[];flash_avg=[];
    
    bar_all=[bar_all;normalized_bar];
    chirp_all=[chirp_all;normalized_chirp];
    color_all=[color_all;normalized_color];
    flash_all=[flash_all;normalized_flash];
    
    
end
% All=[unnormalized_chirp unnormalized_bar unnormalized_color ];% normalizing after groupping Incorrect!!!
% All_Normalized=[(All)'./max(abs(All)')]';

SPK_trg=[trg_data trg_label2 trg_label];% for all cells

 %%

end
