%% 基于四路实测包络的单用户 DOA 估计与周期数测试
% 这个脚本的测试结果也还算可以，但是需要检查参考相位标定和同步采样的问题
%
% 观测模型（单用户、存在相干本地参考/本振）：
%
%   z_m[n] = | h_m + beta_m * exp(j*omega*n) |,
%
% 其中 h_m = alpha*a_m(theta,phi) 是第 m 个阵元的未知复响应。
% 将幅度平方后：
%
%   z_m[n]^2 = |h_m|^2 + |beta_m|^2
%              + 2*Re{h_m*conj(beta_m)*exp(-j*omega*n)}.
%
% 因此，如果实测包络中存在稳定拍频，一个周期内的不同采样点就相当于
% 多个不同参考相位的观测，可以用于相位恢复。重复多个周期主要用于降噪。
%
% 重要限制：
% 1. 四路 ADC 必须同步采样；
% 2. 四路参考链路的相对相位 betaPhaseCalibration 必须已知或完成标定；
% 3. 若包络不是由相干参考形成的稳定拍频，则本脚本输出的 DOA 不可信；
% 4. 实信号无法判断拍频旋转方向，默认会存在 phi 相差约 180 度的镜像解。

clear;
clc;
close all;

%% ======================== 用户配置 ========================

inputFileNames = {
    'adc1_data_amplitude.txt'
    'adc2_data_amplitude.txt'
    'adc3_data_amplitude.txt'
    'adc4_data_amplitude.txt'
};

% 数据分析区间。默认沿用 test2.m 中人工选择的 55%～60%。
analysisPercent = [55, 60];

% 若已知包络拍频，直接填写 Hz；否则保持 NaN，由脚本自动估计。
knownEnvelopeFrequencyHz = NaN;
frequencySearchHz = [0.05, 20];

% 用多少个周期执行主 DOA/GS/GN 测试。
mainCycleCount = 5;

% 用于判断“取几个周期较稳定”的候选周期数。
cycleCountCandidates = [1, 2, 3, 5, 8, 10];
maximumWindowsPerCycleCount = 20;

% 射频和阵列参数。
fc = 3.5e9;
c0 = 3e8;
lambda = c0 / fc;
dx = lambda / 2;
dy = lambda / 2;

% 阵元顺序必须与四个文件一致：左上、右上、左下、右下。
elementPosition = [
    -dx/2,  dy/2
     dx/2,  dy/2
    -dx/2, -dy/2
     dx/2, -dy/2
];

% 四路参考链路的相对复增益。至少必须标定相对相位。
% 未标定时只能先用 ones(4,1) 做可行性测试，不能视为可靠测角结果。
betaPhaseCalibration = ones(4, 1);

% +1：假设参考相位随样本正向旋转；-1：选择共轭镜像解。
% 仅凭实值包络无法判断该符号，应根据收发频率高低或已知来波方向确定。
beatFrequencySign = +1;

% 角度定义与 test2.m 一致：
% theta 为离开阵列法向的角度，phi 为阵列平面内方位角。
thetaGridDeg = 0:1:80;
phiGridDeg = -180:2:180;

% GS/GN 参数。
gsIterationCount = 500;
gnIterationCount = 300;

% 通常本地参考比未知来波强；该假设只用于估计 GS/GN 的参考幅度。
referenceDominates = true;

% 可选真实角度。若未知，保持 NaN。
knownThetaDeg = NaN;
knownPhiDeg = NaN;

%% ======================== 路径和输入检查 ========================

scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    scriptFolder = pwd;
else
    scriptFolder = fileparts(scriptPath);
end

convertedFolder = fullfile(scriptFolder, 'converted_amplitude');
simulationFolder = fullfile(fileparts(scriptFolder), 'simulation');
addpath(simulationFolder);

numberOfElements = length(inputFileNames);
if numberOfElements ~= size(elementPosition, 1)
    error('输入文件数必须等于 elementPosition 的行数。');
end

if length(betaPhaseCalibration) ~= numberOfElements
    error('betaPhaseCalibration 的长度必须等于阵元数。');
end
betaPhaseCalibration = betaPhaseCalibration(:);
betaPhaseCalibration = betaPhaseCalibration ./ ...
    max(abs(betaPhaseCalibration), eps);

timeCell = cell(numberOfElements, 1);
powerDbCell = cell(numberOfElements, 1);

for m = 1:numberOfElements
    inputPath = fullfile(convertedFolder, inputFileNames{m});
    if ~isfile(inputPath)
        error('未找到转换后的幅度文件：%s', inputPath);
    end

    inputTable = readtable(inputPath, ...
        'Delimiter', ',', 'VariableNamingRule', 'preserve');
    if width(inputTable) < 4
        error('%s 至少需要四列。', inputFileNames{m});
    end

    time = inputTable{:, 1};
    powerDb = inputTable{:, 3};
    validIndex = isfinite(time) & isfinite(powerDb);

    timeCell{m} = time(validIndex);
    powerDbCell{m} = powerDb(validIndex);
end

minimumLength = min(cellfun(@length, timeCell));
if minimumLength < 20
    error('有效同步数据点过少。');
end

timeMatrix = zeros(minimumLength, numberOfElements);
powerDbMatrix = zeros(minimumLength, numberOfElements);
for m = 1:numberOfElements
    timeMatrix(:, m) = timeCell{m}(1:minimumLength);
    powerDbMatrix(:, m) = powerDbCell{m}(1:minimumLength);
end

timeMismatchPerSample = max(timeMatrix, [], 2) - ...
    min(timeMatrix, [], 2);
sortedTimeMismatch = sort(timeMismatchPerSample);
timeMismatch95 = sortedTimeMismatch(max(1, ...
    ceil(0.95 * length(sortedTimeMismatch))));
timeMismatch = max(timeMismatchPerSample);
time = median(timeMatrix, 2);
phaseTime = time - time(1);
effectiveSampleRate = (minimumLength - 1) / (time(end) - time(1));

if ~isfinite(effectiveSampleRate) || effectiveSampleRate <= 0
    error('无法从时间列得到有效采样率。');
end

fprintf('\n========== 实测数据基本信息 ==========\n');
fprintf('同步点数                 : %d\n', minimumLength);
fprintf('等效采样率               : %.4f Hz\n', effectiveSampleRate);
fprintf('四路最大时间戳偏差       : %.6g s\n', timeMismatch);
fprintf('四路 95%% 时间戳偏差      : %.6g s\n', timeMismatch95);

if timeMismatch95 > 2 / effectiveSampleRate
    warning(['四路时间戳偏差超过两个等效采样间隔。', ...
        '若 ADC 并非同步采集，跨阵元相位不可用于 DOA。']);
elseif timeMismatch > 2 / effectiveSampleRate
    warning(['存在少量较大的四路时间戳偏差，但 95% 分位偏差正常。', ...
        '脚本使用四路时间戳中位数并保留原始采集顺序。']);
end

% 不能直接使用转换文件中“每路各自最大值归一化”的第四列做跨通道比较。
% 这里由第三列 Power_dB 使用同一个全局参考重新生成线性幅度，保留通道比例。
globalMaximumDb = max(powerDbMatrix(:));
measuredAmplitude = 10.^((powerDbMatrix - globalMaximumDb) / 20);
measuredPower = measuredAmplitude.^2;

%% ======================== 选择分析区间 ========================

if analysisPercent(1) < 0 || analysisPercent(2) > 100 || ...
        analysisPercent(1) >= analysisPercent(2)
    error('analysisPercent 必须满足 0 <= start < end <= 100。');
end

startIndex = floor(minimumLength * analysisPercent(1) / 100) + 1;
endIndex = floor(minimumLength * analysisPercent(2) / 100);
startIndex = max(startIndex, 1);
endIndex = min(endIndex, minimumLength);
analysisIndex = (startIndex:endIndex).';

if length(analysisIndex) < 20
    error('所选分析区间的数据点过少。');
end

%% ======================== 估计共同包络周期 ========================

if isfinite(knownEnvelopeFrequencyHz)
    envelopeFrequencyHz = knownEnvelopeFrequencyHz;
else
    [envelopeFrequencyHz, frequencyFitCost] = ...
        estimateCommonEnvelopeFrequency( ...
        measuredPower(analysisIndex, :), ...
        phaseTime(analysisIndex), effectiveSampleRate, ...
        frequencySearchHz);
    fprintf('自动频率拟合归一化残差   : %.6f\n', frequencyFitCost);
end

samplesPerCycle = effectiveSampleRate / envelopeFrequencyHz;
availableCycleCount = length(analysisIndex) / samplesPerCycle;

fprintf('估计包络频率             : %.6f Hz\n', envelopeFrequencyHz);
fprintf('每周期等效采样点         : %.2f\n', samplesPerCycle);
fprintf('分析区间包含周期数       : %.2f\n', availableCycleCount);

if samplesPerCycle < 6
    warning('每周期少于 6 个采样点，相位扫描采样不足。');
end
if availableCycleCount < 1
    error('分析区间不足一个完整周期。');
end

%% ======================== 正弦模型质量检查 ========================

[~, fullFitInfo] = estimateChannelResponseFromHarmonic( ...
    measuredPower(analysisIndex, :).', phaseTime(analysisIndex), ...
    envelopeFrequencyHz, betaPhaseCalibration, ...
    beatFrequencySign, referenceDominates);

fprintf('\n========== 包络正弦模型检查 ==========\n');
for m = 1:numberOfElements
    fprintf('阵元 %d: 功率正弦拟合 R^2 = %.4f\n', ...
        m, fullFitInfo.rSquared(m));
end
fprintf('四路中位 R^2             : %.4f\n', ...
    median(fullFitInfo.rSquared));

if median(fullFitInfo.rSquared) < 0.5
    warning(['包络的正弦拟合度偏低。当前数据可能不是稳定相干拍频，', ...
        '或对数检波器斜率/采样同步存在问题；DOA 仅供诊断。']);
end

%% ======================== 不同周期数的稳定性扫描 ========================

scanRows = [];
scanWindowResults = struct();

fprintf('\n========== 周期数稳定性测试（直接谐波 ML） ==========\n');

for candidateIndex = 1:length(cycleCountCandidates)
    cycleCount = cycleCountCandidates(candidateIndex);
    windowLength = round(cycleCount * samplesPerCycle);
    possibleWindowCount = floor(length(analysisIndex) / windowLength);
    windowCount = min(possibleWindowCount, maximumWindowsPerCycleCount);

    if windowCount < 1
        continue;
    end

    thetaValues = zeros(windowCount, 1);
    phiValues = zeros(windowCount, 1);
    uValues = zeros(windowCount, 1);
    vValues = zeros(windowCount, 1);
    rSquaredValues = zeros(windowCount, 1);

    for windowIndex = 1:windowCount
        localStart = analysisIndex(1) + (windowIndex - 1) * windowLength;
        localIndex = (localStart:localStart + windowLength - 1).';

        hEstimate = estimateChannelResponseFromHarmonic( ...
            measuredPower(localIndex, :).', phaseTime(localIndex), ...
            envelopeFrequencyHz, betaPhaseCalibration, ...
            beatFrequencySign, referenceDominates);

        [thetaValues(windowIndex), phiValues(windowIndex), ~, ...
            uValues(windowIndex), vValues(windowIndex)] = ...
            estimateDoaFromComplexResponse( ...
            hEstimate, elementPosition, lambda, ...
            thetaGridDeg, phiGridDeg);

        [~, fitInfoWindow] = estimateChannelResponseFromHarmonic( ...
            measuredPower(localIndex, :).', phaseTime(localIndex), ...
            envelopeFrequencyHz, betaPhaseCalibration, ...
            beatFrequencySign, referenceDominates);
        rSquaredValues(windowIndex) = median(fitInfoWindow.rSquared);
    end

    medianTheta = median(thetaValues);
    medianPhi = circularMeanDegrees(phiValues);
    thetaStd = std(thetaValues);
    phiStd = circularStdDegrees(phiValues);
    uStd = std(uValues);
    vStd = std(vValues);
    medianR2 = median(rSquaredValues);

    scanRows = [
        scanRows
        cycleCount, windowCount, medianTheta, medianPhi, ...
        thetaStd, phiStd, uStd, vStd, medianR2
    ]; %#ok<AGROW>

    fieldName = sprintf('cycles_%d', cycleCount);
    scanWindowResults.(fieldName).theta = thetaValues;
    scanWindowResults.(fieldName).phi = phiValues;

    fprintf(['%2d 周期: 窗口数=%2d, theta=%7.2f±%5.2f°, ', ...
        'phi=%8.2f±%5.2f°, R^2中位数=%.3f\n'], ...
        cycleCount, windowCount, medianTheta, thetaStd, ...
        medianPhi, phiStd, medianR2);
end

if isempty(scanRows)
    error('没有候选周期数适合当前分析区间。');
end

cycleScanTable = array2table(scanRows, 'VariableNames', {
    'CycleCount', 'WindowCount', 'MedianTheta_deg', 'MedianPhi_deg', ...
    'StdTheta_deg', 'CircularStdPhi_deg', 'StdU', 'StdV', 'MedianR2'});

cycleResultPath = fullfile(scriptFolder, ...
    'doa_cycle_count_stability.txt');
writetable(cycleScanTable, cycleResultPath, 'Delimiter', ',');

% 无真实角度时，只能按重复窗口的稳定性给建议，不能评价绝对准确度。
stableIndex = find(cycleScanTable.WindowCount >= 2 & ...
    cycleScanTable.StdTheta_deg <= 2 & ...
    cycleScanTable.CircularStdPhi_deg <= 5 & ...
    cycleScanTable.MedianR2 >= 0.5, 1, 'first');

if isempty(stableIndex)
    recommendedCycleCount = NaN;
    fprintf(['未找到同时满足 theta 标准差<=2°、phi 标准差<=5°、', ...
        'R^2>=0.5 的周期数。\n']);
else
    recommendedCycleCount = cycleScanTable.CycleCount(stableIndex);
    fprintf('按重复窗口稳定性，建议至少使用 %d 个周期。\n', ...
        recommendedCycleCount);
end

%% ======================== 主窗口：直接 ML、GS、GN 对比 ========================

actualMainCycleCount = min(mainCycleCount, floor(availableCycleCount));
mainWindowLength = round(actualMainCycleCount * samplesPerCycle);
mainWindowLength = min(mainWindowLength, length(analysisIndex));
mainIndex = analysisIndex(1:mainWindowLength);

[hHarmonic, mainFitInfo] = estimateChannelResponseFromHarmonic( ...
    measuredPower(mainIndex, :).', phaseTime(mainIndex), ...
    envelopeFrequencyHz, betaPhaseCalibration, ...
    beatFrequencySign, referenceDominates);

[thetaHarmonic, phiHarmonic, metricMap, uHarmonic, vHarmonic] = ...
    estimateDoaFromComplexResponse( ...
    hHarmonic, elementPosition, lambda, thetaGridDeg, phiGridDeg);

% 实值包络的共轭镜像候选。
[thetaMirror, phiMirror] = estimateDoaFromComplexResponse( ...
    conj(hHarmonic), elementPosition, lambda, ...
    thetaGridDeg, phiGridDeg);

% 根据功率正弦的直流项和调制度估计参考幅度，供 GS/GN 使用。
betaMagnitude = mainFitInfo.betaMagnitude;
betaEstimate = betaMagnitude .* betaPhaseCalibration;

observationAmplitude = measuredAmplitude(mainIndex, :).';
referenceRotation = exp(1j * beatFrequencySign * 2*pi * ...
    envelopeFrequencyHz * phaseTime(mainIndex).');
referenceMatrix = betaEstimate .* referenceRotation;

% 每个时刻的四个观测依次排列。
zObservation = observationAmplitude(:);
referenceVector = referenceMatrix(:);
observationMatrix = kron(ones(length(mainIndex), 1), ...
    eye(numberOfElements));

fprintf('\n========== 主窗口算法对比（%d 个周期） ==========\n', ...
    actualMainCycleCount);
fprintf('直接谐波 ML : theta=%8.3f°, phi=%8.3f°\n', ...
    thetaHarmonic, phiHarmonic);
fprintf('共轭镜像候选: theta=%8.3f°, phi=%8.3f°\n', ...
    thetaMirror, phiMirror);

gsSucceeded = false;
gnSucceeded = false;
thetaGs = NaN;
phiGs = NaN;
thetaGn = NaN;
phiGn = NaN;

if exist('biased_gs_algorithm', 'file') == 2
    try
        hGs = biased_gs_algorithm( ...
            zObservation, observationMatrix.', ...
            referenceVector, gsIterationCount);
        hGs = normalizeComplexResponse(hGs);
        [thetaGs, phiGs] = estimateDoaFromComplexResponse( ...
            hGs, elementPosition, lambda, thetaGridDeg, phiGridDeg);
        gsSucceeded = true;
        fprintf('biased GS    : theta=%8.3f°, phi=%8.3f°\n', ...
            thetaGs, phiGs);
    catch gsError
        warning('DOA:BiasedGNFailed', ...
                'biased GS 运行失败：%s', gnError.message);
    end
end

if exist('biased_gn_algorithm', 'file') == 2
    try
        hGn = biased_gn_algorithm( ...
            zObservation, observationMatrix, ...
            referenceVector, gnIterationCount);
        hGn = normalizeComplexResponse(hGn);
        [thetaGn, phiGn] = estimateDoaFromComplexResponse( ...
            hGn, elementPosition, lambda, thetaGridDeg, phiGridDeg);
        gnSucceeded = true;
        fprintf('biased GN    : theta=%8.3f°, phi=%8.3f°\n', ...
            thetaGn, phiGn);
    catch gnError
        warning('DOA:BiasedGNFailed', ...
                'biased GN 运行失败：%s', gnError.message);
    end
end

if isfinite(knownThetaDeg) && isfinite(knownPhiDeg)
    fprintf('\n真实角度      : theta=%8.3f°, phi=%8.3f°\n', ...
        knownThetaDeg, knownPhiDeg);
    fprintf('直接 ML 误差 : theta=%+.3f°, phi=%+.3f°\n', ...
        thetaHarmonic - knownThetaDeg, ...
        wrapDegree(phiHarmonic - knownPhiDeg));
end

%% ======================== 保存主结果 ========================

methodNames = {'Harmonic_ML'; 'Conjugate_mirror'; ...
    'Biased_GS'; 'Biased_GN'};
thetaResult = [thetaHarmonic; thetaMirror; thetaGs; thetaGn];
phiResult = [phiHarmonic; phiMirror; phiGs; phiGn];
methodSucceeded = [true; true; gsSucceeded; gnSucceeded];

resultTable = table(methodNames, thetaResult, phiResult, ...
    methodSucceeded, 'VariableNames', ...
    {'Method', 'Theta_deg', 'Phi_deg', 'Succeeded'});
resultPath = fullfile(scriptFolder, 'doa_measured_results.txt');
writetable(resultTable, resultPath, 'Delimiter', ',');

%% ======================== 可视化 ========================

figure('Name', '实测包络及正弦拟合', 'Color', 'w', ...
    'Position', [100, 80, 1200, 760]);
waveformLayout = tiledlayout(2, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

phaseMain = 2*pi * envelopeFrequencyHz * phaseTime(mainIndex);
for m = 1:numberOfElements
    nexttile;
    measuredPowerMain = measuredPower(mainIndex, m);
    fittedPowerMain = mainFitInfo.dc(m) + ...
        mainFitInfo.cosineCoefficient(m) * cos(phaseMain) + ...
        mainFitInfo.sineCoefficient(m) * sin(phaseMain);
    plot(time(mainIndex), measuredPowerMain, '-', 'LineWidth', 0.9);
    hold on;
    plot(time(mainIndex), fittedPowerMain, '-', 'LineWidth', 1.5);
    hold off;
    grid on;
    xlabel('时间 (s)');
    ylabel('相对线性功率');
    title(sprintf('阵元 %d，R^2=%.3f', m, ...
        mainFitInfo.rSquared(m)));
    legend('实测', '正弦拟合', 'Location', 'best');
end
title(waveformLayout, sprintf( ...
    '包络拍频拟合：%.4f Hz，%d 个周期', ...
    envelopeFrequencyHz, actualMainCycleCount));

figure('Name', 'DOA ML 目标函数', 'Color', 'w', ...
    'Position', [180, 100, 1000, 650]);
imagesc(phiGridDeg, thetaGridDeg, ...
    10*log10(metricMap / max(min(metricMap(:)), eps) + eps));
axis xy;
colorbar;
xlabel('方位角 \phi (°)');
ylabel('离开阵列法向角 \theta (°)');
title('直接谐波 ML 归一化残差');
hold on;
plot(phiHarmonic, thetaHarmonic, 'rx', ...
    'MarkerSize', 12, 'LineWidth', 2);
hold off;

figure('Name', '周期数稳定性', 'Color', 'w', ...
    'Position', [220, 130, 950, 520]);
yyaxis left;
plot(cycleScanTable.CycleCount, cycleScanTable.StdTheta_deg, ...
    '-o', 'LineWidth', 1.4);
hold on;
plot(cycleScanTable.CycleCount, ...
    cycleScanTable.CircularStdPhi_deg, ...
    '-s', 'LineWidth', 1.4);
ylabel('跨窗口角度标准差 (°)');
yyaxis right;
plot(cycleScanTable.CycleCount, cycleScanTable.MedianR2, ...
    '-d', 'LineWidth', 1.4);
ylabel('功率正弦拟合 R^2');
xlabel('每个估计窗口使用的周期数');
grid on;
legend('\theta 标准差', '\phi 循环标准差', 'R^2 中位数', ...
    'Location', 'best');
title('周期数与估计稳定性的关系');

fprintf('\n结果已保存：\n');
fprintf('  %s\n', resultPath);
fprintf('  %s\n', cycleResultPath);
if ~isfinite(recommendedCycleCount)
    fprintf(['当前数据没有得到可靠的周期数建议；请优先检查参考相位标定、', ...
        '同步采样和包络正弦拟合度。\n']);
end

%% ======================== 局部函数 ========================

function [frequencyHz, normalizedCost] = ...
        estimateCommonEnvelopeFrequency( ...
        powerData, timeCoordinate, sampleRate, searchHz)

    [sampleCount, channelCount] = size(powerData);
    normalizedData = zeros(size(powerData));
    for channelIndex = 1:channelCount
        y = powerData(:, channelIndex);
        y = detrend(y, 1);
        scale = std(y);
        if scale <= eps
            scale = 1;
        end
        normalizedData(:, channelIndex) = y / scale;
    end

    fftLength = 2^nextpow2(max(4096, 8 * sampleCount));
    window = 0.5 - 0.5*cos(2*pi*(0:sampleCount-1).' / ...
        max(sampleCount - 1, 1));
    spectrum = zeros(fftLength, 1);
    for channelIndex = 1:channelCount
        transform = fft(normalizedData(:, channelIndex) .* ...
            window, fftLength);
        spectrum = spectrum + abs(transform).^2;
    end

    frequencyAxis = (0:fftLength-1).' / fftLength;
    minimumFrequencyCycles = max(searchHz(1) / sampleRate, ...
        1 / (4 * sampleCount));
    maximumFrequencyCycles = min(searchHz(2) / sampleRate, 0.45);
    searchIndex = find(frequencyAxis >= minimumFrequencyCycles & ...
        frequencyAxis <= maximumFrequencyCycles);

    if isempty(searchIndex)
        error('frequencySearchHz 与数据采样率不匹配。');
    end

    [~, localPeakIndex] = max(spectrum(searchIndex));
    initialFrequencyHz = frequencyAxis(searchIndex(localPeakIndex)) * ...
        sampleRate;

    objective = @(frequency) commonSinusoidCost( ...
        frequency, normalizedData, timeCoordinate, ...
        searchHz(1), searchHz(2));
    options = optimset('Display', 'off', 'TolX', 1e-10, ...
        'MaxIter', 300, 'MaxFunEvals', 1000);
    frequencyHz = fminsearch( ...
        objective, initialFrequencyHz, options);
    frequencyHz = min(max(frequencyHz, searchHz(1)), searchHz(2));
    normalizedCost = objective(frequencyHz);
end

function cost = commonSinusoidCost( ...
        frequency, data, timeCoordinate, ...
        minimumFrequency, maximumFrequency)

    if frequency < minimumFrequency || frequency > maximumFrequency
        distance = max(minimumFrequency - frequency, 0) + ...
            max(frequency - maximumFrequency, 0);
        cost = 1e6 * (1 + distance^2);
        return;
    end

    timeCoordinate = timeCoordinate(:);
    localTime = timeCoordinate - timeCoordinate(1);
    phase = 2*pi*frequency*localTime;
    designMatrix = [ones(size(localTime)), ...
        localTime / max(localTime(end), eps), ...
        cos(phase), sin(phase)];

    residualEnergy = 0;
    totalEnergy = 0;
    for channelIndex = 1:size(data, 2)
        y = data(:, channelIndex);
        coefficient = designMatrix \ y;
        residual = y - designMatrix * coefficient;
        residualEnergy = residualEnergy + sum(residual.^2);
        totalEnergy = totalEnergy + sum(y.^2);
    end
    cost = residualEnergy / max(totalEnergy, eps);
end

function [hEstimate, info] = estimateChannelResponseFromHarmonic( ...
        powerMatrix, timeCoordinate, frequency, betaCalibration, ...
        beatSign, referenceDominates)

    % powerMatrix: M x N
    timeCoordinate = timeCoordinate(:);
    phase = 2*pi*frequency*timeCoordinate;
    designMatrix = [ones(length(timeCoordinate), 1), ...
        cos(phase), sin(phase)];

    elementCount = size(powerMatrix, 1);
    dc = zeros(elementCount, 1);
    cosineCoefficient = zeros(elementCount, 1);
    sineCoefficient = zeros(elementCount, 1);
    rSquared = zeros(elementCount, 1);
    betaMagnitude = zeros(elementCount, 1);
    hMagnitude = zeros(elementCount, 1);
    crossTerm = zeros(elementCount, 1);

    for m = 1:elementCount
        y = powerMatrix(m, :).';
        coefficient = designMatrix \ y;
        fitted = designMatrix * coefficient;

        dc(m) = coefficient(1);
        cosineCoefficient(m) = coefficient(2);
        sineCoefficient(m) = coefficient(3);
        residualEnergy = sum((y - fitted).^2);
        totalEnergy = sum((y - mean(y)).^2);
        rSquared(m) = 1 - residualEnergy / max(totalEnergy, eps);

        crossTerm(m) = (cosineCoefficient(m) + ...
            1j * beatSign * sineCoefficient(m)) / 2;

        modulation = hypot(cosineCoefficient(m), ...
            sineCoefficient(m));
        discriminant = max(dc(m)^2 - modulation^2, 0);
        largerPower = max((dc(m) + sqrt(discriminant)) / 2, eps);
        smallerPower = max((dc(m) - sqrt(discriminant)) / 2, eps);

        if referenceDominates
            betaMagnitude(m) = sqrt(largerPower);
            hMagnitude(m) = sqrt(smallerPower);
        else
            betaMagnitude(m) = sqrt(smallerPower);
            hMagnitude(m) = sqrt(largerPower);
        end
    end

    betaEstimate = betaMagnitude .* betaCalibration(:);
    betaEstimate(abs(betaEstimate) < eps) = eps;
    hEstimate = crossTerm ./ conj(betaEstimate);
    hEstimate = normalizeComplexResponse(hEstimate);

    info.dc = dc;
    info.cosineCoefficient = cosineCoefficient;
    info.sineCoefficient = sineCoefficient;
    info.rSquared = rSquared;
    info.betaMagnitude = betaMagnitude;
    info.hMagnitude = hMagnitude;
    info.crossTerm = crossTerm;
end

function hNormalized = normalizeComplexResponse(h)
    h = h(:);
    hNormalized = zeros(size(h));
    valid = abs(h) > 1e-12 & isfinite(h);
    hNormalized(valid) = h(valid) ./ abs(h(valid));
end

function [thetaEstimate, phiEstimate, metricMap, ...
        uEstimate, vEstimate] = estimateDoaFromComplexResponse( ...
        hEstimate, elementPosition, lambda, thetaGrid, phiGrid)

    hEstimate = normalizeComplexResponse(hEstimate);
    if norm(hEstimate) <= eps
        error('复阵列响应接近零，无法估计 DOA。');
    end

    [phiMesh, thetaMesh] = meshgrid(phiGrid, thetaGrid);
    uGrid = sind(thetaMesh(:)) .* cosd(phiMesh(:));
    vGrid = sind(thetaMesh(:)) .* sind(phiMesh(:));
    waveNumber = 2*pi / lambda;

    steeringMatrix = exp(-1j * waveNumber * ( ...
        elementPosition(:, 1) * uGrid.' + ...
        elementPosition(:, 2) * vGrid.'));

    coherence = abs(steeringMatrix' * hEstimate).^2 ./ ...
        max(sum(abs(steeringMatrix).^2, 1).' * norm(hEstimate)^2, eps);
    metricVector = 1 - coherence;
    metricMap = reshape(metricVector, length(thetaGrid), ...
        length(phiGrid));

    [~, bestIndex] = min(metricVector);
    thetaEstimate = thetaMesh(bestIndex);
    phiEstimate = phiMesh(bestIndex);
    uEstimate = uGrid(bestIndex);
    vEstimate = vGrid(bestIndex);
end

function meanDegree = circularMeanDegrees(angleDegree)
    meanComplex = mean(exp(1j * deg2rad(angleDegree)));
    meanDegree = rad2deg(angle(meanComplex));
end

function stdDegree = circularStdDegrees(angleDegree)
    resultantLength = abs(mean(exp(1j * deg2rad(angleDegree))));
    resultantLength = min(max(resultantLength, eps), 1);
    stdDegree = rad2deg(sqrt(-2*log(resultantLength)));
end

function wrappedDegree = wrapDegree(angleDegree)
    wrappedDegree = mod(angleDegree + 180, 360) - 180;
end
