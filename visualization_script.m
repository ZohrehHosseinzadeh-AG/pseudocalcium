% This code plots the clusterd RGC data 

%Note: This code requiers:
% * RGC data from Baden et al (https://doi.org/10.5061/dryad.d9v38). 
% * SPCA toolbox to run SPCA algorithm of Zou et. al[1] for
% computing sparse principal components (DOI: 10.18637/jss.v084.i10).  
% * jbfill function from John A. Bockstege (November 2006) for visualzing data.


%    [1] H. Zou, T. Hastie, and R. Tibshirani. Sparse Principal Component
%    Analysis. J. Computational and Graphical Stat. 15(2):265-286, 2006.

%% Read Baden data and extract features

clear
close all

load('D:\Hame2\Data\Baden\BadenEtAl_RGCs_2016_v1.mat')
cluster_idx(cluster_idx==-1)=76;% set noisy clusters as last one
for cl=1:length(unique(cluster_idx))
    
    m_chirp_baden(cl,:)= mean(chirp_avg(:,cluster_idx==cl),2);
    m_bar_baden(cl,:)= mean(bar_tc(:,cluster_idx==cl),2);
    m_color_baden(cl,:)= mean(color_avg(:,cluster_idx==cl),2);
      
end

all_baden=[m_chirp_baden m_bar_baden m_color_baden];
all_baden=[m_chirp_baden m_bar_baden m_color_baden];
[B1, SD1, L, D] = spca(chirp_avg', [] , 20, inf,-10);%10-Non-Zero Time-bins p<<n will result in inf

%Getting 20 features for each cluster projected on principle component
%directions
P_chirp_baden = B1'* m_chirp_baden';%Matrix containing 20 features per cluster(For Chirp Stimulus)

%Time Course response to Moving bar (Got from SVD)
[B2, SD2] = spca(bar_tc', [] , 8, inf,-5);%5-Non-Zero Time-bins p<<n will result in inf
P_bar_baden = B2'* m_bar_baden';%Matrix containing 8 features per cluster(For Moving-bar Stimulus)

%Other 4 features should be taken from it's derivatcie
% first derivative , Step size assumed to be 1
[B3, SD3] = spca(diff(bar_tc)', [] , 4, inf,-6);%6-Non-Zero Time-bins p<<n will result in inf
P_der_bar_baden = B3'* m_bar_baden(:,1:end-1)';%Matrix containing 8 features per cluster(For Moving-bar Stimulus)

[B4, SD4] = spca(color_avg', [] , 6, inf,-10);%10-Non-Zero Time-bins p<<n will result in inf
P_color_baden = B4'*m_color_baden';%Matrix containing 6 features per cluster(For Moving-bar Stimulus)

A_C_F=[P_chirp_baden;P_bar_baden;P_der_bar_baden;P_color_baden];

%% Read our data and average them



load('.\Data\data_list.mat')% run Metadata_list 
[SPK_org,PSTH,var_idx,chirp_org, bar_org, color_org,flash_org,data,nr_name_org,Bad_cell_name,e_STAs_org]=averaging_func(ds_list);% averaging data

%% find bad cell indx
%
e_STAs=e_STAs_org;
nr_name=nr_name_org;
bar_all=bar_org;
chirp_all=chirp_org;
color_all=color_org;
flash_all=flash_org;
SPK_all=SPK_org;
var_idx(isnan(var_idx))=1;
nodatacelles=find(isnan(var_idx));
all=[chirp_all bar_all color_all];

for ii=1:length(nodatacelles)
    SPK_all(end+1,3)=nodatacelles(ii);
end

rmv_cells_l=(var_idx>=1);% cells with variance ratio higehr than 1
bad_cellall=sum(isnan(all),2)>1;
bad_cellall_bar=find(sum(isnan(bar_all),2)>1);



%% compute distance matrix
SPK_all(:,4)=SPK_all(:,3);
lbl=[]
[a f]=unique(SPK_all(:,3));
for ii=2:length(f)
    lbl=[lbl; ii-1*ones(f(ii)-f(ii-1),1)];
end
lbl(end+1:size(SPK_all,1))=(ii);
SPK_all(:,3)=lbl;


SPK_all(:,4)=SPK_all(:,3);
lbl=[];
[a f]=unique(SPK_all(:,3));
for ii=2:length(f)
    lbl=[lbl; ii-1*ones(f(ii)-f(ii-1),1)];
end
lbl(end+1:size(SPK_all,1))=(ii);
SPK_all(:,3)=lbl;

normalized_chirp=chirp_all;
normalized_bar=bar_all;
normalized_color=color_all;
normalized_flash=flash_all;

%% Project data onto extracted principle components
P_Cell_Chirp=B1'*(normalized_chirp');
P_Cell_MV = B2'*(normalized_bar');
P_Cell_MV_Der = B3'*(diff(normalized_bar'));
P_Cell_Color =B4'*(normalized_color');
T_C_F=[P_Cell_Chirp;P_Cell_MV;P_Cell_MV_Der;P_Cell_Color];


for k=1:size(A_C_F,2)
    for j=1:size(T_C_F,2)
        
        [Corr_Coef,P_val]=corrcoef(T_C_F(:,j),A_C_F(:,k));% computes correlation and p-values
        
        Dist(k,j) = norm(T_C_F(:,j)- A_C_F(:,k));% distance matirx
        Dist_man(k,j) = sqrt(sum((T_C_F(:,j)-A_C_F(:,k)).^2)/numel(A_C_F(:,k)));
        
        CR(k,j)= Corr_Coef(1,2);
        PV(k,j)=P_val(1,2);
    end
end
Dist=Dist_man;

[cr_val,idx]=min(Dist);
%subplot(313)

for ii=1:size(PV,2)
    pvs(ii)=PV(idx(ii),ii);
end
%scatter(idx,pvs)
%% significance of peak and trough of STA

for clnam=1:size(nr_name,1)
    nr_name{clnam,2}=idx(clnam);
    e_STAs(clnam).ch=idx(clnam);
    
end
e_stas=[]
for stas=1:length(e_STAs)
    if ~isempty(e_STAs(stas).esta)
        e_stas(stas,:)=e_STAs(stas).esta;
    else
        e_stas(stas,:)=zeros(1,50);
    end
end


n_e_stas=(e_stas-mean(e_stas,2))./std(e_stas')';
n_e_stas(isnan(n_e_stas))=0;
for n_test=1:length(e_STAs)
    STA=e_stas(n_test,:)
    baseline_sta = mean(STA(length(STA)/2+10:end));
    std_baselines = std(STA(length(STA)/2+10:end));
    tKerLen=25
    alphaville = 1-nthroot(.95,tKerLen);
    
    
    alphaville = alphaville/2;
    peak_sta = min(STA(1:length(STA)/2));
    if sum (STA~=0)
        H_peak(n_test) = ztest(peak_sta,baseline_sta,std_baselines,alphaville);
        
        trough_sta = max(STA(1:length(STA)/2));
        H_trough(n_test) = ztest(trough_sta,baseline_sta,std_baselines,alphaville);
        %with_sta(n_test) = kstest(n_e_stas(n_test,:));
    else
        H_peak(n_test)=0;
        H_trough(n_test)=0;
        with_sta(n_test) = 0;
    end
end
with_sta=H_trough|H_peak;% cells with significant STA

sta_t=linspace(-1,1,50);
sum(with_sta)

sta_f=figure
subplot(122)
imagesc(sta_t,[],n_e_stas(with_sta,:))
title('Significant STAs')
ylabel('Number of cells')
xlabel('Time (s)')

colorbar
subplot(121)
imagesc(sta_t,[],n_e_stas(~with_sta,:))
title('Non-significant STAs')
xlabel('Time (s)')
colorbar
print(sta_f,'bad_good_stas','-dpng')
%% Cells to be used for final analysis
Condition=[with_sta'&~bad_cellall&~rmv_cells_l];
%Condition=[~bad_cellall&~rmv_cells_l];

sum(Condition)
n_e_stas=n_e_stas(Condition,:);
idx=idx(Condition);
Dist=Dist(:,Condition);
normalized_chirp=normalized_chirp(Condition,:);
normalized_bar=normalized_bar(Condition,:);
normalized_color=normalized_color(Condition,:);
normalized_flash=normalized_flash(Condition,:);
%% histogram of clusters
cluster_idx_baden=cluster_idx;
cluster_idx_baden(cluster_idx_baden>38)=[];

dist_cluster=figure
idx_no76=idx;
idx_no76(idx_no76>38)=[];
histogram(cluster_idx_baden,[0:39])
hold on
histogram(idx_no76,[0:39])
xlabel('Cluster number')
xticks([1:38])

ylabel('Number of cells')
legend('Reference data','Our data')
%print(dist_cluster,'dist_cluster','-dpng')
xlim([1,39])
grid on




%% comparing distance matrix

[cluster_idx_my,sidx]=sort(idx);
nDist=(Dist./max(Dist));
figure
imagesc(Dist')
title('Distance matrix before clustering')
ylabel('Cell #')
xlabel('Cluster #')
colorbar

figure
imagesc(Dist(:,sidx)')
title('Distance matrix after clustering')
ylabel('Cell #')
xlabel('Cluster #')
colorbar

%plot(crs(sidx))

figure
imagesc(CR')
title('Corelation coeff before clustering')
ylabel('Cell #')
xlabel('Cluster #')
colorbar
figure
imagesc(CR(:,sidx)')
title('Corelation coeff after clustering')
ylabel('Cell #')
xlabel('Cluster #')
colorbar
%% sort data by cluster number
%chirp_avg_clusterd=zeros(size(A_C_F,2),size(chirp_avg,2));
chirp_avg_clusterd=normalized_chirp(sidx,:);
m_bar_clusterd=normalized_bar(sidx,:);
color_clusterd=normalized_color(sidx,:);
flash_clusterd=normalized_flash(sidx,:);
All_clusterd=all(sidx,:);
hold on
histogram(idx,76)
subplot(312)
scatter(cluster_idx_my,cr_val(sidx))
title('corr coef value for each cell in each cluster')




%% plot mean clusters
colors=jet(49);
colors=colors(linspace(49,1,49),:)
colors(50:76,:)=.5*ones(76-49,3)
aa=(linspace(76,1,76))
bb=(linspace(1,76,76))

hamed_hist=histc(idx,[0:76])
hamed_hist(1)=[]
ratios=100*hamed_hist/sum(hamed_hist(1:49))
for i=aa
    ylbl{i}=[num2str(aa(i)) ,'    (' ,num2str(ceil(ratios(aa(i)))) '%)  ']
end

fig_cluster=figure('units','normalized','outerposition',[0 0 .4 3]);
%subplot(121)
noise=zeros(size(normalized_chirp,1),16)

hold on

patch([0,50,50,0],[76-50,76-50,77,77],[0 0 0],'Marker','.','EdgeColor',[1 1 1])


alpha .05
for cl=1:76
    strt=76+1;
    % if sum((idx==cl))>0
    cluster_number(cl)=sum((idx==cl))
    
    m_sta(cl,:)= nanmean(e_stas(idx==cl,:),1);
    
    minn=min([median(normalized_chirp(idx==cl,:),1) median(normalized_bar(idx==cl,:),1)  median(noise(idx==cl,:),1) median(normalized_color(idx==cl,:),1)]);
    
    m_all_hamed(cl,:)=[mean(normalized_chirp(idx==cl,:),1) mean(normalized_bar(idx==cl,:),1)...
        nanmean(normalized_color(idx==cl,:),1)];
    
    
    
    [ph,msg]=jbfill(linspace(0,48,length(m_all_hamed(cl,:))),...
        .8*m_all_hamed(cl,:)-cl+strt, +.8*(min(m_all_hamed(cl,:)))+...
        zeros(1,length(m_all_hamed(cl,:)))-cl+strt,colors(cl,:),[0 0 0],0,1);
    set(gca,'YTick',[1:76],'YTickLabel',ylbl,'FontSize',12,'Layer','top','XTick',[0 32 36 42 48],...
        'XTickLabel',{'0','32','36','42','48'})
    hold on
    
    %    end
end



line([2,2],[76-50 80],'LineStyle','-.','Color',[.5 .5 .5])
line([5,5],[76-50 80],'LineStyle','-.','Color',[.5 .5 .5])
line([8,8],[76-50 80],'LineStyle','-.','Color',[.5 .5 .5])
line([10,10],[76-50 80],'LineStyle','-.','Color',[.5 .5 .5])
line([18,18],[76-50 80],'LineStyle','-.','Color',[.5 .5 .5])
line([20,20],[76-50 80],'LineStyle','-.','Color',[.5 .5 .5])
line([28,28],[76-50 80],'LineStyle','-.','Color',[.5 .5 .5])
line([32,32],[76-50 80],'LineStyle','-.','Color',[.5 .5 .5])
line([36,36],[76-50 80],'LineStyle','-.','Color',[.5 .5 .5])
line([39,39],[76-50 80],'LineStyle','-.','Color',[.5 .5 .5])
line([39+3,39+3],[76-50 80],'LineStyle','-.','Color',[.5 .5 .5])
line([42+3,42+3],[76-50 80],'LineStyle','-.','Color',[.5 .5 .5])
yl= ylabel('Cluster Number')
xlabel('time (s)')
set(gca,'FontSize',10)
ylim([27.5,77])
xlim([0,48])
%    print(fig_cluster,'fig_cluster_pca_significant_STA_less_data','-dpng')
%    print(fig_cluster,'fig_cluster_pca_significant_STA_less_data','-deps')
%   saveas(fig_cluster,'fig_cluster_pca_significant_STA_less_data')
%% plot stas

close all
clear m_n_e_stas
f=figure('units','normalized','outerposition',[0 0 .2 1.5]);
adds=[1:76];
patch([-1,2,2,-1],[76-50,76-50,-77,-77],[.95 .95 .95],'Marker','.','EdgeColor',[1 1 1])

hold on

%n_e_stas=(e_stas-mean(e_stas,2))./std(e_stas')';
n_e_stas(isnan(n_e_stas))=0;
n_e_stas2=n_e_stas/2
x=linspace(-1,1,length(n_e_stas2(1,:)))

for cl=1:76
    
    m_n_e_stas(cl,:)=nanmedian(n_e_stas2(idx==cl,:),1)-adds(cl);
    plot(x,m_n_e_stas(cl,:),'Color',colors(cl,:),'LineWidth',3)
    line([-1,1],[mean(m_n_e_stas(cl,:)),mean(m_n_e_stas(cl,:))],'Color','black')
    
    hold on
    
end
box off
set(gca,'ytick',[])
ylim([-50,0])
xlim([-1,.5])
a = get(gca,'XTickLabel');
xlabel('time (s)')

set(gca,'XTickLabel',a,'FontName','Times','fontsize',7)
line([0,0],[-70,0],'LineStyle','-.','Color',[.5 .5 .5])
% print(f,'fig_cluster_esta_significant_less_data','-dpng')
% print(f,'fig_cluster_esta_significant_less_data','-deps')
% saveas(f,'fig_cluster_esta_significant_less_data')
