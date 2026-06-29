%% 这个脚本尝试对adc绘制的电压信号进行处理，能显示波形形状，目前结果比较理想
clear;
clc;
close all;

%% ================= 文件设置 =================
fileNames = {
    'adc1_data.txt'
    'adc2_data.txt'
    'adc3_data.txt'
    'adc4_data.txt'
};

figureTitles = {
    '文件1'
    '文件2'
    '文件3'
    '文件4'
};

%% ================= 局部数据范围设置 =================
% 按各文件自身的数据点总数截取相同比例
startPercent = 12;      % 起始百分比
endPercent   = 18;      % 结束百分比

if startPercent < 0 || endPercent > 100 || startPercent >= endPercent
    error('百分比范围应满足：0 <= startPercent < endPercent <= 100。');
end

%% ================= 数据存储 =================
numberOfFiles = length(fileNames);

timeData = cell(numberOfFiles, 1);
voltageData = cell(numberOfFiles, 1);

%% ================= 读取四个文件 =================
for k = 1:numberOfFiles

    fileName = fileNames{k};

    if ~isfile(fileName)
        warning('未找到文件：%s', fileName);
        continue;
    end

    % 自动识别文件格式
    opts = detectImportOptions(fileName, ...
        'Delimiter', ',', ...
        'VariableNamingRule', 'preserve');

    dataTable = readtable(fileName, opts);

    % 第一列为时间，第二列为电压
    time = dataTable{:, 1};
    voltage = dataTable{:, 2};

    % 去除无效数据
    validIndex = isfinite(time) & isfinite(voltage);

    time = time(validIndex);
    voltage = voltage(validIndex);

    % 不建议按照时间重新排序
    % 保持文件原始采集顺序，避免相同时间点之间的顺序被改变
    timeData{k} = time;
    voltageData{k} = voltage;

    fprintf('文件%d：%s\n', k, fileName);
    fprintf('有效数据点总数：%d\n\n', length(voltage));

end

%% ================= 第一张图：四组完整数据 =================
figure('Position', [100, 80, 1200, 760]);

fullLayout = tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for k = 1:numberOfFiles

    nexttile;

    time = timeData{k};
    voltage = voltageData{k};

    if isempty(time)
        axis off;

        text(0.5, 0.5, ...
            ['无有效数据：', fileNames{k}], ...
            'HorizontalAlignment', 'center');

        continue;
    end

    plot(time, voltage, '-', ...
        'LineWidth', 1.0);

    grid on;
    box on;

    xlabel('时间 (s)');
    ylabel('电压 (V)');
    title(figureTitles{k});

    if min(time) < max(time)
        xlim([min(time), max(time)]);
    end

end

title(fullLayout, '四组完整电压数据');

%% ================= 第二张图：四组局部数据 =================
figure('Position', [130, 100, 1200, 760]);

localLayout = tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

for k = 1:numberOfFiles

    nexttile;

    time = timeData{k};
    voltage = voltageData{k};

    if isempty(time)
        axis off;

        text(0.5, 0.5, ...
            ['无有效数据：', fileNames{k}], ...
            'HorizontalAlignment', 'center');

        continue;
    end

    %% 根据百分比计算索引
    totalPointNumber = length(voltage);

    startIndex = floor(totalPointNumber * startPercent / 100) + 1;
    endIndex = floor(totalPointNumber * endPercent / 100);

    % 防止越界
    startIndex = max(startIndex, 1);
    endIndex = min(endIndex, totalPointNumber);

    if endIndex < startIndex
        axis off;

        text(0.5, 0.5, ...
            '所选区间内数据点不足', ...
            'HorizontalAlignment', 'center');

        continue;
    end

    %% 提取局部数据
    localTime = time(startIndex:endIndex);
    localVoltage = voltage(startIndex:endIndex);

    %% 绘图
    plot(localTime, localVoltage, '-', ...
        'LineWidth', 1.0);

    grid on;
    box on;

    xlabel('时间 (s)');
    ylabel('电压 (V)');

    title(sprintf('%s：%.2f%%～%.2f%%', ...
        figureTitles{k}, startPercent, endPercent));

    if min(localTime) < max(localTime)
        xlim([min(localTime), max(localTime)]);
    end

    %% 输出区间信息
    fprintf('文件%d局部区间：\n', k);
    fprintf('总数据点数：%d\n', totalPointNumber);
    fprintf('数据序号：%d ～ %d\n', startIndex, endIndex);
    fprintf('局部数据点数：%d\n\n', length(localVoltage));

end

title(localLayout, sprintf( ...
    '四组局部电压数据：%.2f%%～%.2f%%', ...
    startPercent, endPercent));