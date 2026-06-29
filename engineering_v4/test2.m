%% 这个脚本尝试对接收到的电压数据信号进行分析，通过角度域的建模完成doa估计，假设单参考相位
% 结果还算理想，之后可以拿这套算法多测几组实验
clear;
clc;
close all;

%% =========================================================
%  1. 四个阵元文件
%
%  顺序：
%  1：左上
%  2：右上
%  3：左下
%  4：右下
%% =========================================================
fileNames = {
    'adc1_data.txt'
    'adc2_data.txt'
    'adc3_data.txt'
    'adc4_data.txt'
};

elementNames = {
    '左上'
    '右上'
    '左下'
    '右下'
};
figureTitles = {
    '左上'
    '右上'
    '左下'
    '右下'
};


numberOfElements = 4;

%% 先展示完整数据图形
% ================= 数据存储 =================
numberOfFiles = length(fileNames);

timeData = cell(numberOfFiles, 1);
voltageData = cell(numberOfFiles, 1);
% ================= 读取四个文件 =================
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

%% =========================================================
%  2. 选择需要分析的数据区间
%% =========================================================
startPercent = 55;
endPercent   = 60;

if startPercent < 0 || endPercent > 100 || ...
        startPercent >= endPercent
    error('百分比范围应满足 0 <= startPercent < endPercent <= 100。');
end

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

%% =========================================================
%  3. 射频和阵列参数
%% =========================================================
fc = 3.5e9;
c0 = 3e8;

lambda = c0/fc;
k0 = 2*pi/lambda;

% 阵元间距，根据实际尺寸修改
dx = lambda/2;
dy = lambda/2;

% 阵列位于x-y平面，法向为z轴
%
% 顺序：左上、右上、左下、右下
elementPosition = [
    -dx/2,  dy/2;
     dx/2,  dy/2;
    -dx/2, -dy/2;
     dx/2, -dy/2
];

%% =========================================================
%  4. 对数检波器的电压—功率标定参数
%
%  假设：
%  Vout = voltageOffset + voltageSlope * PdB
%
%  voltageSlope单位为V/dB。
%
%  下面的0.025 V/dB只是示例，必须根据实际检波器修改。
%% =========================================================
voltageOffset = 0;
voltageSlope  = 0.025;

%% =========================================================
%  5. 四路本振复增益标定值
%
%  betaCalibration(m)包含：
%  1. 功分器各输出端的幅度差；
%  2. 馈线长度引起的相位差；
%  3. 合路器和接收链路的固定幅相差。
%
%  未标定时先全部设为1。
%% =========================================================
betaCalibration = ones(4,1);

% 标定后可以写成如下形式：
%
% betaCalibration = [
%     1.00 * exp(1j*deg2rad(0));
%     0.97 * exp(1j*deg2rad(8));
%     1.03 * exp(1j*deg2rad(-6));
%     0.99 * exp(1j*deg2rad(4))
% ];

%% =========================================================
%  6. 读取并处理四路数据
%% =========================================================
timeData = cell(numberOfElements,1);
voltageData = cell(numberOfElements,1);
powerData = cell(numberOfElements,1);

measuredPower = zeros(numberOfElements,1);

for m = 1:numberOfElements

    fileName = fileNames{m};

    if ~isfile(fileName)
        error('未找到文件：%s', fileName);
    end

    opts = detectImportOptions(fileName, ...
        'Delimiter', ',', ...
        'VariableNamingRule', 'preserve');

    dataTable = readtable(fileName, opts);

    time = dataTable{:,1};
    voltage = dataTable{:,2};

    validIndex = isfinite(time) & isfinite(voltage);

    time = time(validIndex);
    voltage = voltage(validIndex);

    totalNumber = length(voltage);

    startIndex = floor(totalNumber*startPercent/100) + 1;
    endIndex   = floor(totalNumber*endPercent/100);

    startIndex = max(startIndex,1);
    endIndex   = min(endIndex,totalNumber);

    if endIndex < startIndex
        error('%s所选区间中没有足够的数据。', elementNames{m});
    end

    time = time(startIndex:endIndex);
    voltage = voltage(startIndex:endIndex);

    %% 对数检波电压恢复为相对线性功率
    %
    % PdB = (Vout - offset)/slope
    %
    % Plinear ∝ 10^(PdB/10)
    %
    relativePower = 10.^((voltage-voltageOffset) ...
                         /(10*voltageSlope));

    timeData{m} = time;
    voltageData{m} = voltage;
    powerData{m} = relativePower;

    % 对所选区间取中位数，抑制异常点
    measuredPower(m) = median(relativePower);

    fprintf('%s：总数据点数 = %d，分析区间点数 = %d\n', ...
        elementNames{m}, totalNumber, length(voltage));
end

%% =========================================================
%  7. 功率归一化
%
%  公共包络功率和检波器的公共比例因子可以通过归一化消除。
%% =========================================================
if sum(measuredPower) <= 0
    error('测得功率异常。');
end

measuredPowerNormalized = measuredPower/sum(measuredPower);

fprintf('\n归一化实测功率：\n');

for m = 1:numberOfElements
    fprintf('%s：%.6f\n', ...
        elementNames{m}, measuredPowerNormalized(m));
end

%% =========================================================
%  8. 非线性参数估计
%
%  优化参数：
%  p(1) = u
%  p(2) = v
%  p(3) = log(rho)，其中rho为来波相对于本振的幅度
%  p(4) = 来波公共初相位psi
%
%  方向余弦：
%  u = sin(theta)cos(phi)
%  v = sin(theta)sin(phi)
%% =========================================================
objectiveFunction = @(p) directionObjective( ...
    p, ...
    measuredPowerNormalized, ...
    elementPosition, ...
    betaCalibration, ...
    k0);

% 使用多个随机初值，降低陷入局部极小值的概率
numberOfStarts = 100;

bestObjective = inf;
bestParameter = [];

rng(1);

options = optimset( ...
    'Display', 'off', ...
    'MaxIter', 2000, ...
    'MaxFunEvals', 5000, ...
    'TolX', 1e-10, ...
    'TolFun', 1e-12);

for n = 1:numberOfStarts

    % u、v初值位于单位圆附近
    initialU = 2*rand-1;
    initialV = 2*rand-1;

    if initialU^2 + initialV^2 > 1
        scale = sqrt(initialU^2 + initialV^2);
        initialU = 0.9*initialU/scale;
        initialV = 0.9*initialV/scale;
    end

    % 来波相对本振幅度初值
    initialLogRho = log(0.1 + 2*rand);

    % 公共相位初值
    initialPsi = -pi + 2*pi*rand;

    initialParameter = [
        initialU
        initialV
        initialLogRho
        initialPsi
    ];

    [estimatedParameter, objectiveValue] = fminsearch( ...
        objectiveFunction, ...
        initialParameter, ...
        options);

    if objectiveValue < bestObjective
        bestObjective = objectiveValue;
        bestParameter = estimatedParameter;
    end
end

%% =========================================================
%  9. 参数恢复
%% =========================================================
estimatedU = bestParameter(1);
estimatedV = bestParameter(2);

estimatedRho = exp(bestParameter(3));
estimatedPsi = wrapToPiLocal(bestParameter(4));

directionRadius = sqrt(estimatedU^2 + estimatedV^2);

% 数值误差保护
directionRadius = min(max(directionRadius,0),1);

estimatedTheta = asind(directionRadius);
estimatedPhi = atan2d(estimatedV,estimatedU);

fprintf('\n================ 估计结果 ================\n');
fprintf('u = %.6f\n', estimatedU);
fprintf('v = %.6f\n', estimatedV);
fprintf('离开阵列法向的角度 theta = %.2f°\n', estimatedTheta);
fprintf('阵列平面方位角 phi       = %.2f°\n', estimatedPhi);
fprintf('来波/本振幅度比 rho      = %.4f\n', estimatedRho);
fprintf('来波公共相位 psi         = %.2f°\n', rad2deg(estimatedPsi));
fprintf('归一化拟合误差           = %.6e\n', bestObjective);

%% =========================================================
%  10. 根据估计参数计算拟合功率
%% =========================================================
phaseVector = -k0*( ...
    elementPosition(:,1)*estimatedU + ...
    elementPosition(:,2)*estimatedV);

steeringVector = exp(1j*phaseVector);

alphaEstimated = estimatedRho*exp(1j*estimatedPsi);

fittedPower = abs( ...
    alphaEstimated*steeringVector + betaCalibration ...
    ).^2;

fittedPowerNormalized = fittedPower/sum(fittedPower);

% %% =========================================================
% %  11. 四路局部波形
% %% =========================================================
% figure('Position',[100,80,1200,760]);
% 
% waveformLayout = tiledlayout(2,2, ...
%     'TileSpacing','compact', ...
%     'Padding','compact');
% 
% for m = 1:numberOfElements
% 
%     nexttile;
% 
%     plot(timeData{m},voltageData{m},'-', ...
%         'LineWidth',1.0);
% 
%     grid on;
%     box on;
% 
%     xlabel('时间 (s)');
%     ylabel('检波电压 (V)');
%     title(elementNames{m});
% end
% 
% title(waveformLayout,sprintf( ...
%     '四阵元局部数据：%.2f%%～%.2f%%', ...
%     startPercent,endPercent));

%% =========================================================
%  12. 实测功率和拟合功率比较
%% =========================================================
figure('Position',[250,150,850,480]);

barData = [
    measuredPowerNormalized, ...
    fittedPowerNormalized
];

bar(barData);

grid on;
box on;

xticks(1:4);
xticklabels(elementNames);

ylabel('归一化功率');
legend('实测功率','模型拟合功率', ...
    'Location','best');

title(sprintf( ...
    '固定本振测角拟合：\\theta=%.2f°，\\phi=%.2f°', ...
    estimatedTheta,estimatedPhi));

%% =========================================================
%  13. 绘制目标函数的二维局部角度谱
%
%  固定已经估计出的rho和psi，只扫描u、v对应的角度。
%% =========================================================
thetaGrid = 0:1:80;
phiGrid = -180:2:180;

angleCost = zeros(length(thetaGrid),length(phiGrid));

for it = 1:length(thetaGrid)

    thetaTrial = thetaGrid(it);

    for ip = 1:length(phiGrid)

        phiTrial = phiGrid(ip);

        uTrial = sind(thetaTrial)*cosd(phiTrial);
        vTrial = sind(thetaTrial)*sind(phiTrial);

        phaseTrial = -k0*( ...
            elementPosition(:,1)*uTrial + ...
            elementPosition(:,2)*vTrial);

        steeringTrial = exp(1j*phaseTrial);

        predictedPower = abs( ...
            alphaEstimated*steeringTrial + ...
            betaCalibration ...
            ).^2;

        predictedPower = predictedPower/sum(predictedPower);

        angleCost(it,ip) = sum( ...
            (measuredPowerNormalized-predictedPower).^2);
    end
end

figure('Position',[180,100,1000,650]);

imagesc(phiGrid,thetaGrid, ...
    10*log10(angleCost/min(angleCost(:)) + eps));

axis xy;
colorbar;

xlabel('方位角 \phi (°)');
ylabel('离开阵列法向的角度 \theta (°)');
title('固定本振功率匹配代价函数');

hold on;

plot(estimatedPhi,estimatedTheta,'rx', ...
    'MarkerSize',12, ...
    'LineWidth',2);

hold off;

%% =========================================================
%  局部函数
%% =========================================================
function cost = directionObjective( ...
    parameter, measuredPower, ...
    elementPosition, betaCalibration, k0)

    u = parameter(1);
    v = parameter(2);

    rho = exp(parameter(3));
    psi = parameter(4);

    %% 方向余弦必须满足u^2+v^2<=1
    radiusSquared = u^2 + v^2;

    if radiusSquared > 1
        directionPenalty = 1e3*(radiusSquared-1)^2;
    else
        directionPenalty = 0;
    end

    phaseVector = -k0*( ...
        elementPosition(:,1)*u + ...
        elementPosition(:,2)*v);

    steeringVector = exp(1j*phaseVector);

    alpha = rho*exp(1j*psi);

    predictedPower = abs( ...
        alpha*steeringVector + ...
        betaCalibration ...
        ).^2;

    powerSum = sum(predictedPower);

    if ~isfinite(powerSum) || powerSum <= eps
        cost = 1e10;
        return;
    end

    predictedPower = predictedPower/powerSum;

    fittingError = sum( ...
        (measuredPower-predictedPower).^2);

    cost = fittingError + directionPenalty;
end

function wrappedPhase = wrapToPiLocal(phase)
    wrappedPhase = mod(phase+pi,2*pi)-pi;
end