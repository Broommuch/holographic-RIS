%% 这一版是文件的初版处理，只能实现图形的显示功能

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

% 每张子图的标题，可根据实验条件修改
figureTitles = {
    '文件1'
    '文件2'
    '文件3'
    '文件4'
};

%% ================= 创建画布 =================
figure('Position', [200, 100, 1100, 750]);

t = tiledlayout(2, 2, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

%% ================= 依次读取并绘图 =================
for k = 1:length(fileNames)

    fileName = fileNames{k};

    % 检查文件是否存在
    if ~isfile(fileName)
        warning('未找到文件：%s', fileName);

        nexttile;
        axis off;
        text(0.5, 0.5, ...
            ['未找到文件：', fileName], ...
            'HorizontalAlignment', 'center');

        continue;
    end

    %% 读取数据
    opts = detectImportOptions(fileName, ...
        'Delimiter', ',', ...
        'VariableNamingRule', 'preserve');

    dataTable = readtable(fileName, opts);

    % 第一列为时间，第二列为电压
    time = dataTable{:, 1};
    voltage = dataTable{:, 2};

    %% 去除无效数据
    validIndex = isfinite(time) & isfinite(voltage);

    time = time(validIndex);
    voltage = voltage(validIndex);

    %% 按时间排序
    [time, sortIndex] = sort(time);
    voltage = voltage(sortIndex);

    %% 绘制子图
    nexttile;

    plot(time, voltage, '-', ...
        'LineWidth', 1.2);

    grid on;
    box on;

    xlabel('时间 (s)');
    ylabel('电压 (V)');
    title(figureTitles{k});

    % 数据不为空时设置横坐标范围
    if ~isempty(time) && min(time) < max(time)
        xlim([min(time), max(time)]);
    end

end

%% ================= 整体标题 =================
title(t, '四组电压数据随时间变化曲线');