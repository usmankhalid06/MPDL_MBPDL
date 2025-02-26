clear;
close all; 
clc;

%% Initializing parameters
var_arr = [1 0.05];
x1 = 21; 
x2 = 21; 
nV = x1; 
nSRCS = 6; 
N = 250; 
nIter = 30; %% iterations for dictionary
RT = 1;
K = 9;

%% Constructing spatial and temporal sources
warning('off')
onset_intervals = [50, 45, 50, 30, 40, 60];
durations = [15, 18, 21, 12, 20, 30];
TC = zeros(N, length(onset_intervals)); % Preallocate matrix for efficiency

for i = 1:length(onset_intervals)
    onsets = 0:onset_intervals(i):N-20;
    TC(:, i) = generate_TC(onsets, durations(i), N);
end


% hrf = spm_hrf(RT);
% R  = [TC(1,2) zeros(1,length(hrf)-1)];
for i = 1:size(TC, 2)
%     TC(:, i) = zscore(toeplitz(TC(:, i), R)*hrf); % Apply HRF and toeplitz structure if needed
    TC(:, i) = zscore(TC(:, i) + 0.005 * wfbm(1, N)'); % Add fractional Brownian motion if needed
end

Dp = dctbases(N, N); % DCT basis dictionary


numSlices = 6; % Number of slices
tmpSM = zeros(x1, x2, numSlices); % Initialize all slices at once

% Define regions and Gaussian sigma for each slice
regions = {
    [2, 8, 2, 8, 5]    % [startRow, endRow, startCol, endCol, sigma]
    [2, 9, 11, 19, 5]
    [8, 13, 2, 11, 7]
    [8, 15, 11, 19, 9]
    [12, 19, 2, 11, 10]
    [13, 19, 11, 19, 6]
};

% Apply Gaussian activation regions to each slice
for i = 1:numSlices
    region = regions{i};
    rows = region(2) - region(1) + 1;
    cols = region(4) - region(3) + 1;
    sigma = region(5);
    % Create a Gaussian mask
    G = gaussian2D(rows, cols, sigma);
    % Place the Gaussian mask into the slice
    tmpSM(region(1):region(2), region(3):region(4), i) = G;
    SM(i, :) = reshape(tmpSM(:, :, i), 1, x1*x2); % Flatten and assign to SM
    SM(i,:) = zscore(SM(i,:));
end


%% Generating data
Y = (TC+sqrt(var_arr(1))*randn(N,nSRCS))*(SM+sqrt(var_arr(2))*randn(nSRCS,x1*x2)); 
Y= Y-repmat(mean(Y),size(Y,1),1);


%% ODL  
% for ODL to work you need to download SPAMS toolbox (thats what we used in paper)
% however, when sharing the code on github we made a custom mexLasso as part of ODL function
tStart=tic;
spa = 15; %15
bs = 3; nIter2 = 30; %size(Y,2)/bs;
[D{1},X{1},E1,C1]= my_ODL(Y, K, nIter,spa,bs,nIter2,TC,SM);
tEnd(1,1) = toc(tStart);
% figure; plot(E1)
% figure; plot(C1)

%% ACSD
spa = 100;
tStart=tic;
[D{2},X{2},E2,C2]= my_ACSD(Y,Dp(:,1:K),spa,nIter,TC,SM); 
tEnd(1,1) = toc(tStart);
% figure; plot(E2)
% figure; plot(C2)

%% ACSDBE
spa = 100;
lam = 30;
tStart=tic;
[D{3},X{3},E3,C3]= my_ACSDBE(Y,Dp(:,1:150),K,spa,lam,nIter,TC,SM); 
tEnd(1,1) = toc(tStart);
% figure; plot(E3)
% figure; plot(C3)

%% MPDL
PSF_sim3 = spcol(augknt(linspace(0,1,35),4),4,linspace(0,1,21))'; 
PSF_sim4 = kron(PSF_sim3, PSF_sim3); 
PSF_sim4 = abs(PSF_sim4); 
for jjj=1:size(PSF_sim4,1)
    PSF_sim4(jjj,:) =(PSF_sim4(jjj,:) - min(PSF_sim4(jjj,:))) / ( max(PSF_sim4(jjj,:)) - min(PSF_sim4(jjj,:)) );
end
spa = 100;
lam = [30,90];
tStart=tic;
[D{4},X{4},E4,C4]= my_MPDL(Y,Dp(:,1:150),PSF_sim4,K,spa,lam(1),lam(2),nIter,TC,SM);
tEnd(1,1) = toc(tStart);
% figure; plot(E4)
% figure; plot(C4)


%% MBPDL
tStart=tic;
PSF_sim3 = spcol(augknt(linspace(0,1,35),4),4,linspace(0,1,21))'; 
PSF_sim4 = kron(PSF_sim3, PSF_sim3); 
PSF_sim4 = abs(PSF_sim4); 
for jjj=1:size(PSF_sim4,1)
    PSF_sim4(jjj,:) =(PSF_sim4(jjj,:) - min(PSF_sim4(jjj,:))) / ( max(PSF_sim4(jjj,:)) - min(PSF_sim4(jjj,:)) );
end
SM2 = abs(SM); 
for jjj=1:size(SM2,1)
    SM2(jjj,:) =(SM2(jjj,:) - min(SM2(jjj,:))) / ( max(SM2(jjj,:)) - min(SM2(jjj,:)) );
end
f = 4; % asuming only first 3 sources are known
l1 = 5;
l2 = 18;
l3 = 30;
l4 = 90;
TC2 =  TC*diag(1./sqrt(sum(TC.*TC))); 
spa = 100;
[D{5},X{5}, E5,C5]= my_MBPDL(Y,[TC2(:,1:4) Dp(:,1:150)],[SM2(1:4,:); PSF_sim4],K,spa,l1,l2,l3,l4,nIter,f,TC,SM); % and 2 %=
tEnd(1,2) = toc(tStart)
% figure; plot(E5)
% figure; plot(C5)


%% Analysis
nA=6;
sD{1} = TC;
sX{1} = SM;
for jj =1:nA-1
    [sD{jj+1},sX{jj+1},ind]=sort_TSandSM_spatial(TC,SM,D{jj},X{jj},nSRCS);
    for ii =1:nSRCS
        TCcorr(jj+1,ii) =abs(corr(TC(:,ii),D{jj}(:,ind(ii))));
        SMcorr(jj+1,ii) =abs(corr(SM(ii,:)',X{jj}(ind(ii),:)'));
    end
end



TP = zeros(nA,nSRCS); FP = zeros(nA,nSRCS); FN = zeros(nA,nSRCS); F_score = 0; thr = 0.05;
for i =1:nA-1
    for jjj=1:nSRCS
        SM_gw4(jjj,:) =SM(jjj,:)/norm(SM(jjj,:)); 
        [~, indd(jjj)]  = max(abs(corr(abs(SM_gw4(jjj,:)'),abs(X{i}'))));
        XX{i}(indd(jjj),:) = X{i}(indd(jjj),:)/norm(X{i}(indd(jjj),:));
        TP(i,jjj) = TP(i,jjj) +sum(abs(SM_gw4(jjj,:))>=thr & abs(XX{i}(indd(jjj),:))>=thr);
        FP(i,jjj) = FP(i,jjj) +sum(abs(SM_gw4(jjj,:))<=thr & abs(XX{i}(indd(jjj),:))>=thr);
        FN(i,jjj) = FN(i,jjj) +sum(abs(SM_gw4(jjj,:))>=thr & abs(XX{i}(indd(jjj),:))<=thr); 
        Fscore(i,jjj) = (2*sum(TP(i,jjj)))/(2*sum(TP(i,jjj))+sum(FP(i,:))+sum(FN(i,jjj)));
    end
%     F_score(i) = (2*sum(TP(i,:)))/(2*sum(TP(i,:))+sum(FP(i,:))+sum(FN(i,:)));
end
Fscore = [zeros(1,nSRCS); Fscore];


TC_corr = mean(TCcorr');
SM_corr = mean(SMcorr');
F_score = mean(Fscore');
 
f = figure; f.Position = [170 120 1200 600]; my_subplots_fig4_new(nA, nSRCS,nV,nV,TCcorr(:,:,1),Fscore,sD,sX,nSRCS); 


%% Printing results
% Print the header
fprintf('Following are the results for ODL, ACSD, ACSDBE, MPDL, MBPDL:\n');
fprintf('--------------------------------------------------------------------------------------------------\n');

% Table header
fprintf('%-60s', 'Metric');
fprintf('%8s   %8s   %8s   %8s   %8s\n', 'ODL', 'ACSD', 'ACSDBE', 'MPDL', 'MBPDL');

% Correlation between source TCs and recovered TCs
fprintf('%-60s', 'Correlation between source TCs and recovered TCs:');
fprintf('  %8.4f  %8.4f  %8.4f    %8.4f  %8.4f\n', TC_corr(2:end));

% Correlation between source SMs and recovered SMs
fprintf('%-60s', 'Correlation between source SMs and recovered SMs:');
fprintf('  %8.4f  %8.4f  %8.4f    %8.4f  %8.4f\n', SM_corr(2:end));

% Fscores for recovered SMs
fprintf('%-60s', 'Fscores for recovered SMs:');
fprintf('  %8.4f  %8.4f  %8.4f    %8.4f  %8.4f\n', F_score(2:end));

fprintf('--------------------------------------------------------------------------------------------------\n');


