# WaveletFilter_PCA_ERPs
Combine wavelet filter and PCA and use them on ERPs
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
