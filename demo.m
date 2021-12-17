% Demo script for running RIRreconstADMM.m
%
% Required files: RIRreconstADMM.m, RIR_FDTD.mat
% 
% Please refer to the following paper for details:
%
% Kohei Yatabe and Akiko Sugahara
% "Simulation of room impulse response using frequency-domain FEM and
%  convex-optimization-based post-processing"
% Applied Acoustics (2022)

%% parameter setting

cutoffL = 100; % [Hz] 
cutoffH = 800; % [Hz]

strengthL1 = 10;


%% get impulse response

load RIR_FDTD

t = (0:length(h)-1)'/fs;
f = (0:length(h)-1)'*fs/length(h);

[~,maxIdx] = max(abs(h));
preIdx = maxIdx - 1;


%% create data

trueH = fft(h);
lenH = length(trueH);
trueIdx = 1:lenH;

availableIdx = trueIdx(cutoffL<f & f<cutoffH);
data = trueH(availableIdx);


%% zero-filling method

missingX = zeros(lenH,1);
missingX(availableIdx) = data;

zeroFill = ifft(missingX,'symmetric');


%% proposed method

tic
estimated = RIRreconstADMM(data,lenH,availableIdx,preIdx,strengthL1);
toc


%% plot result

figure
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

nexttile
plot(t,h,'LineWidth',1)
xlim([0 0.05])
ylim([-1 1])
xlabel('Time [s]','FontSize',11)
title('True RIR','FontSize',12)

nexttile
plot(t,zeroFill,'LineWidth',1)
xlim([0 0.05])
ylim([-1 1])
xlabel('Time [s]','FontSize',11)
title('Zero-filling method','FontSize',12)

nexttile
plot(t,estimated,'LineWidth',1)
xlim([0 0.05])
ylim([-1 1])
xlabel('Time [s]','FontSize',11)
title('Proposed method','FontSize',12)