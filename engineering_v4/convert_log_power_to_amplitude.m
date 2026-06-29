%% 这个脚本能将对数检波器输出电压转换为线性幅度包络
% 已测试完成
%
% 输入文件的第一列为时间(s)，第二列为对数检波器输出电压(V)。
% 标定模型：
%
%   Vout = voltageOffset + voltageSlope * PdB
%
% 因而：
%
%   PdB              = (Vout - voltageOffset) / voltageSlope
%   linearPower      ∝ 10^(PdB/10)
%   linearAmplitude  ∝ 10^(PdB/20)
%
% 注意：功率检波器只能恢复非负的幅度包络，不能恢复射频信号的正负、
% 载波波形或相位。未做绝对功率标定时，只应使用相对幅度。

clear;
clc;
close all;

%% ======================== 用户配置 ========================

% 脚本会处理 experiment 文件夹中的所有匹配文件。
inputFilePattern = 'adc*_data.txt';

% 对数检波器标定参数。0.025 V/dB 目前沿用 test2.m 中的示例值，
% 必须根据所用检波器的数据手册或实测标定结果修改。
voltageOffset = 0;       % V；PdB = 0 dB 时的检波输出电压
voltageSlope  = 0.025;   % V/dB；下降型检波器应填负值

% false：只生成相对幅度（推荐用于尚未完成绝对标定的数据）。
% true ：同时把 PdB 解释为 dBm，生成负载上的绝对 Vrms 和峰值。
enableAbsoluteConversion = false;
referenceImpedance = 50; % ohm，仅绝对转换使用

% 转换后的 TXT/CSV 文件保存到此子文件夹。
outputFolderName = 'converted_amplitude';

%% ======================== 参数检查 ========================

if ~isfinite(voltageSlope) || voltageSlope == 0
    error('voltageSlope 必须是非零有限数值。');
end

if enableAbsoluteConversion && ...
        (~isfinite(referenceImpedance) || referenceImpedance <= 0)
    error('referenceImpedance 必须为正数。');
end

scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    scriptFolder = pwd;
else
    scriptFolder = fileparts(scriptPath);
end

inputFiles = dir(fullfile(scriptFolder, inputFilePattern));
if isempty(inputFiles)
    error('未找到输入文件：%s', ...
        fullfile(scriptFolder, inputFilePattern));
end

outputFolder = fullfile(scriptFolder, outputFolderName);
if ~isfolder(outputFolder)
    mkdir(outputFolder);
end

%% ======================== 批量转换 ========================

figure('Name', '对数检波电压转换为线性幅度', ...
    'Color', 'w', 'Position', [100, 100, 1100, 700]);
plotLayout = tiledlayout(length(inputFiles), 1, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for fileIndex = 1:length(inputFiles)
    inputPath = fullfile(inputFiles(fileIndex).folder, ...
        inputFiles(fileIndex).name);

    importOptions = detectImportOptions(inputPath, ...
        'Delimiter', ',', 'VariableNamingRule', 'preserve');
    inputTable = readtable(inputPath, importOptions);

    if width(inputTable) < 2
        error('文件至少需要两列（时间和检波电压）：%s', inputPath);
    end

    time = inputTable{:, 1};
    detectorVoltage = inputTable{:, 2};

    validIndex = isfinite(time) & isfinite(detectorVoltage);
    time = time(validIndex);
    detectorVoltage = detectorVoltage(validIndex);

    if isempty(time)
        warning('文件中没有有效数据，已跳过：%s', inputPath);
        continue;
    end

    % 检波电压 -> 对数功率
    powerDb = (detectorVoltage - voltageOffset) ./ voltageSlope;

    % 先减去全文件最大值再做指数运算，使结果最大值为1并避免溢出。
    relativePowerDb = powerDb - max(powerDb);
    relativeAmplitude = 10.^(relativePowerDb ./ 20);
    relativePower = relativeAmplitude.^2;

    resultTable = table(time, detectorVoltage, powerDb, ...
        relativeAmplitude, relativePower, ...
        'VariableNames', {'Time_s', 'DetectorVoltage_V', 'Power_dB', ...
        'RelativeAmplitude', 'RelativePower'});

    if enableAbsoluteConversion
        % 仅当 powerDb 经过标定、单位确实为 dBm 时成立。
        powerW = 1e-3 .* 10.^(powerDb ./ 10);
        amplitudeVrms = sqrt(powerW .* referenceImpedance);
        amplitudeVpeak = sqrt(2) .* amplitudeVrms;

        resultTable.Power_W = powerW;
        resultTable.Amplitude_Vrms = amplitudeVrms;
        resultTable.Amplitude_Vpeak = amplitudeVpeak;
    end

    [~, baseName, ~] = fileparts(inputFiles(fileIndex).name);
    outputPath = fullfile(outputFolder, ...
        [baseName, '_amplitude.csv']);
    writetable(resultTable, outputPath);

    nexttile;
    plot(time, relativeAmplitude, 'LineWidth', 1.0);
    grid on;
    ylabel('相对幅度');
    title(strrep(inputFiles(fileIndex).name, '_', '\_'));
    if fileIndex == length(inputFiles)
        xlabel('时间 (s)');
    end

    fprintf('已转换：%s\n', inputFiles(fileIndex).name);
    fprintf('  有效点数：%d\n', length(time));
    fprintf('  输出文件：%s\n\n', outputPath);
end

title(plotLayout, '线性幅度包络（每路最大值归一化为 1）');

fprintf('转换完成。\n');
if ~enableAbsoluteConversion
    fprintf(['当前输出为相对幅度；取得检波器的 dBm 标定参数后，', ...
        '可开启绝对幅度转换。\n']);
end
