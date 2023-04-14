% Here, the Wavelet filter and Principal component analysis (PCA) have been done with simulated ERP data (20 subjects and two stimulus types)
% The components before and after back-projected are also shown here
%
% The wavelet filter is constructed following the guidance in the book 'Advanced Signal Processing on ERPs', Chapter 2, by Cong Fengyu 
% Cong Fengyu wrote the PCA rotation codes.
%
% Before running this script 'm_waveletFilter_PCA.m', the EEGLAB is necessary, and the 'topoplot' function was called from it
%
%
% simudata.mat; % include A,S,E, correponding to the temporal matrix, spatial matrix, and matrix of background errors
% the simulated data includes 6components, in which N4 changes amplitude between stimuli/conditions, but other components do not 
% compName = {'N1','N2','N4','P1','P2','P3'};
% NumSub = 20;% number of subjects
% NumCond = 2; % number of stimuli/conditions

% temporal information includes fs, NumSamps, timeStart, timeEnd
% fs : sampling frequency of the simulated data, here is 500 Hz
% NumSamps : number of temporal samples of the simulated data, here is 600
% timeStart : here is 200, ERP started from 200ms before stimulus onset
% timeEnd : here is 1000, ERP ended 1000ms after stimulus onset
% 
% spatial information includes interested channel index
% to each component (IntereChan), channel locations, 
% IntereChan = [6,62,24,22,39,58]; % channel indexes, corresponding to the components' peak values of position, components are 'N1','N2','N4','P1','P2','P3' 
% chanlocs = readlocs('P3_N4.elp'); % includes 65 channels
% NumChan = 65;% number of channels 

% X_ERP = A*S; % ERP data reconstruction by temporal and spatial matrixes.
% 
% LevelNoise=1; % the level of noise
% X = X_ERP+E*0.1*LevelNoise; % add noise (in matrix E) to the data
% X = X-ones(NumSamps,1)*mean(X(1:100,:));%  baseline correction
% 
% wavelet filter on ERPs
% wname = 'rbio6.8'; 
% lv = 9;% defined as lv, where fs = 2^lv
% kp = [3 4 5 6];
% for colIndex = 1:size(X_temporal,2)
%     X_temporal(:,colIndex) =
%     f_filterWavelet(X_temporal(:,colIndex),lv,wname,kp); 
% end
% 
% %% PCA 
% [coeff,score,lambda] = princomp(X_temporal');
% 
% K = 8;% the number of components retained
%    
% V = coeff(:,1:K);
% Z = score(:,1:K)';
% 
% %% rotation of retained components
% [V_rotated,transf_V] = rotatefactors(V,'method','promax','power',3,'maxit',500);
% Y = Z'*transf_V; % using structure matrix but not pattern matrix
% Y = Y';
% T = V*inv(transf_V'); 
% Y = reshape(Y,K,NumChan,NumSub,NumCond);

% by Yang Tiantian, written with MATLAB 2016b, in Jyvaskyla University
%%
clear all
close all;
clc
tic
%%
fs = 500;% Hz
NumSamps = 600;
timeStart = -200;%ms
timeEnd = 1000;%ms
tIndex = linspace(timeStart,timeEnd,NumSamps);
% compTimePoint = [30,33,26,26,39,25];%point
IntereChan = [6,62,24,22,39,58];
compName = {'N1','N2','N4','P1','P2','P3'};
chanlocs = readlocs('P3_N4.elp');
NumChan = 65;
NumSub = 20;
NumCond = 2;

%% Data
load simudata; % include A,S,E
X_ERP = A*S;
%% add noise to the data
LevelNoise=1;
X = X_ERP+E*0.1*LevelNoise;
X = X-ones(NumSamps,1)*mean(X(1:100,:));%  baseline correction

X_temporal = reshape(X,NumSamps,NumChan*NumSub*NumCond);
%% wavelet filter
wname = 'rbio6.8';
lv = 9;% defined as lv, where fs = 2^lv
kp = [3 4 5 6];

for colIndex = 1:size(X_temporal,2)
    X_temporal(:,colIndex) = f_filterWavelet(X_temporal(:,colIndex),lv,wname,kp);

%     X_temporal(:,colIndex) = X_temporal(:,colIndex) - mean(X_temporal(1:100,colIndex));
end

%% PCA

[coeff,score,lambda] = princomp(X_temporal');
%% 
K = 8;% the number of components
   
V = coeff(:,1:K);
Z = score(:,1:K)';

%% rotation of retained components
[V_rotated,transf_V] = rotatefactors(V,'method','promax','power',3,'maxit',500);
Y = Z'*transf_V; % using structure matrix but not pattern matrix
Y = Y';
%%
T = V*inv(transf_V'); 
Y = reshape(Y,K,NumChan,NumSub,NumCond);
%%
route = {['PCA_Results_github\retain',int2str(K),' components\']};
mkdir(route{1});
save([route{1} 'PCA'],'T','Y','lambda','K','-v7.3');

%%  Check each component's condition
stiName = {'sti-1','sti-2'};
NumFig = ceil(K/3);
NumPlot = 3;
count = 0;
for fige = 1:NumFig       
    h(fige) = figure;
    set(gcf,'outerposition',get(0,'screensize'));
    for numplot = 1:NumPlot
        count = count+1;
        if count > K
            break;
        end
        subplot(NumPlot,5,1+(numplot-1)*5,'align');

        set(gca,'fontsize',14);
        plot(tIndex,T(:,count),'k','linewidth',3);
        grid on;
        set(gca,'ydir','reverse');
        xlim([timeStart timeEnd]);
        ylabel(['Comp #',int2str(count)]);       
        xlabel(['Time (ms)']);

        y_comp = squeeze(Y(count,:,:,:));
        y_mean = squeeze(mean(y_comp,2));
        maxV = max(abs(y_mean(:)));
        % %
        for sti = 1:NumCond

            stiCorrPlot = (numplot-1)*5+sti+3;

            subplot(NumPlot,5,stiCorrPlot);
            set(gca,'fontsize',14);
            sti_y = squeeze(y_comp(:,:,sti));
            rho = corr(sti_y);
            imagesc(rho);
            set(gca,'clim',[-1 1]);
            title(stiName{sti});
            if sti == 2
                xlabel('Subject #');
                hs = subplot(NumPlot,5,stiCorrPlot);
                hs_position = get(hs,'position');
                hc = colorbar('peer',gca,'fontsize',14);
                set(get(hc,'ylabel'),'string','Correlation coefficient','fontsize',14);
                set(hs,'pos',hs_position);
            end

            stiTopoPlot = stiCorrPlot-2;

            subplot(NumPlot,5,stiTopoPlot);
            set(gca,'fontsize',14);
            topo = mean(sti_y,2);
            topoplot(topo,chanlocs);
            caxis([-maxV maxV]);
            title(stiName{sti});
            if sti == 1 % fix the colorbar
                hs =  subplot(NumPlot,5,stiTopoPlot);
                hs_position = get(hs,'position');
                hc = colorbar('peer',gca,'fontsize',14);
                set(get(hc,'ylabel'),'string','Amplitude','fontsize',14);
                set(hs,'pos',hs_position);
            end               
        end

    end      
end

%% back projection 
% the first 6 components are the interested ones, so we back=projected them
% and further analyze them separately
% here we are interested in the component #1 which is corresponding to N4,
% and component #6 and #7 are also similar to the N4, thus we can
% back-projected component #1, or combined these three components and
% back-projected them.

IntereCompIndex = [1,6,7];
Y = reshape(Y,K,NumChan*NumSub*NumCond);
compN4 = T(:,IntereCompIndex(1))*Y(IntereCompIndex(1),:);
combinedCompN4 = T(:,IntereCompIndex)*Y(IntereCompIndex,:);

%% check the back-projected component (N4)
h(fige+1) = figure;
set(gcf,'outerposition',get(0,'screensize'));
for ii = 1:2
    backProjectedData = [];
    waveform_channel = [];
    waveform = [];
    if ii == 1
        backProjectedData = reshape(compN4,NumSamps,NumChan,NumSub,NumCond);
        compName = 'N4 (one component back-projected)';
    else
        backProjectedData = reshape(combinedCompN4,NumSamps,NumChan,NumSub,NumCond);
        compName = 'N4 (three component back-projected)';
    end
    waveform_channel = squeeze(backProjectedData(:,IntereChan(3),:,:));%IntereChan(3)is the channel corresponds to N4 
    waveform = squeeze(mean(waveform_channel,2));
    
    % plot waveform
    subplot(2,3,1+(ii-1)*3,'align');
    set(gca,'fontsize',14);
    plot(tIndex,waveform(:,1),'k','linewidth',3);
    grid on;
    hold on;
    plot(tIndex,waveform(:,2),'color',[0.5,0.5,0.5],'linewidth',3);
    set(gca,'ydir','reverse');
    xlim([timeStart timeEnd]);          
    xlabel(['Time (ms)']);
    ylabel('Amplitude (\muV)');
    legend(stiName{1},stiName{2});
    title(compName);

    % plot topography, prepare data
    [peak,col] = max(max(abs(waveform)));
    [peakPoint,col] = find(peak == abs(waveform));
    Topo_on_Peak = squeeze(backProjectedData(peakPoint,:,:,:));
    Topo_meansub = squeeze(mean(Topo_on_Peak,2));
    
    maxV = max(max(abs(Topo_meansub)));
    % %
    for sti = 1:NumCond

        stiTopoPlot = (ii-1)*3+sti+1;
        subplot(2,3,stiTopoPlot,'align');
        set(gca,'fontsize',14);
               
        topoplot(Topo_meansub(:,sti),chanlocs);
        caxis([-maxV maxV]);
        title(stiName{sti});
        if sti == 2 % fix the colorbar
            hs =  subplot(2,3,stiTopoPlot);
            hs_position = get(hs,'position');
            hc = colorbar('peer',gca,'fontsize',14);
            set(get(hc,'ylabel'),'string','Amplitude','fontsize',14);
            set(hs,'pos',hs_position);
        end               
    end
end
%%
toc