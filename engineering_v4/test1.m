%% 这个脚本尝试对接收到的电压数据信号进行分析，通过角度域的建模完成doa估计，因为假设了多参考相位，目前有点小问题
clear;
clc;
close all;

%% =========================================================
%  1. 文件设置
%  每个文件依次包含参考相位：
%  0°、90°、180°、270°
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

numberOfElements = 4;

%% =========================================================
%  2. 射频和阵列参数
%% =========================================================
fc = 3.2e9;                 % 载波频率
c0 = 3e8;
lambda = c0 / fc;
k0 = 2*pi / lambda;

% 阵元间距
dx = lambda/2;
dy = lambda/2;

% 四个阵元坐标
% 顺序：左上、右上、左下、右下
elementPosition = [
    -dx/2,  dy/2;
     dx/2,  dy/2;
    -dx/2, -dy/2;
     dx/2, -dy/2
];

%% =========================================================
%  3. 对数检波器标定参数
%
%  假设：
%  Vout = voltageOffset + voltageSlope * PdB
%
%  voltageSlope单位为V/dB。
%
%  下面只是示例值，必须根据芯片数据手册或实测标定修改。
%% =========================================================
voltageOffset = 0;          % V
voltageSlope = 0.025;       % V/dB，例如25 mV/dB

%% =========================================================
%  4. 数据区间设置
%
%  可以只使用文件中某一百分比区间，然后再将其分成四段。
%% =========================================================
startPercent = 10;
endPercent   = 90;

if startPercent < 0 || endPercent > 100 || ...
        startPercent >= endPercent
    error('百分比范围设置错误。');
end

%% =========================================================
%  5. 每个相位区间去掉边缘过渡部分
%
%  例如discardRatio=0.1表示每段前后各丢掉10%。
%% =========================================================
discardRatio = 0.10;

%% =========================================================
%  6. 读取四个阵元的数据
%% =========================================================
timeData = cell(numberOfElements,1);
voltageData = cell(numberOfElements,1);
powerData = cell(numberOfElements,1);

for m = 1:numberOfElements

    if ~isfile(fileNames{m})
        error('未找到文件：%s', fileNames{m});
    end

    opts = detectImportOptions(fileNames{m}, ...
        'Delimiter', ',', ...
        'VariableNamingRule', 'preserve');

    dataTable = readtable(fileNames{m}, opts);

    time = dataTable{:,1};
    voltage = dataTable{:,2};

    validIndex = isfinite(time) & isfinite(voltage);

    time = time(validIndex);
    voltage = voltage(validIndex);

    % 保持原始采样顺序
    totalNumber = length(voltage);

    startIndex = floor(totalNumber * startPercent/100) + 1;
    endIndex   = floor(totalNumber * endPercent/100);

    startIndex = max(startIndex,1);
    endIndex   = min(endIndex,totalNumber);

    time = time(startIndex:endIndex);
    voltage = voltage(startIndex:endIndex);

    %% 将检波电压恢复成相对线性功率
    %
    % PdB = (Vout - offset)/slope
    % Plinear ∝ 10^(PdB/10)
    %
    powerLinear = 10.^((voltage - voltageOffset) ...
                       /(10*voltageSlope));

    timeData{m} = time;
    voltageData{m} = voltage;
    powerData{m} = powerLinear;

    fprintf('%s阵元：有效数据点数 = %d\n', ...
        elementNames{m}, length(voltage));
end

%% =========================================================
%  7. 显示四个阵元的原始检波电压
%% =========================================================
figure('Position',[100,80,1200,760]);

voltageLayout = tiledlayout(2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

for m = 1:numberOfElements

    nexttile;

    plot(timeData{m}, voltageData{m}, '-', ...
        'LineWidth',1.0);

    grid on;
    box on;

    xlabel('时间 (s)');
    ylabel('检波电压 (V)');
    title(elementNames{m});
end

title(voltageLayout,'四阵元对数检波输出');

%% =========================================================
%  8. 将每个文件平均分成四个参考相位区间
%
%  顺序：
%  q=1：0°
%  q=2：90°
%  q=3：180°
%  q=4：270°
%% =========================================================
meanPower = zeros(numberOfElements,4);

for m = 1:numberOfElements

    powerLinear = powerData{m};
    numberOfSamples = length(powerLinear);

    segmentLength = floor(numberOfSamples/4);

    if segmentLength < 10
        error('%s阵元数据太少，无法分成四个区间。', ...
            elementNames{m});
    end

    for q = 1:4

        segmentStart = (q-1)*segmentLength + 1;

        if q < 4
            segmentEnd = q*segmentLength;
        else
            segmentEnd = numberOfSamples;
        end

        segmentPower = powerLinear(segmentStart:segmentEnd);

        % 去除相位切换时的前后过渡区
        numberInSegment = length(segmentPower);
        discardNumber = floor(discardRatio * numberInSegment);

        validStart = discardNumber + 1;
        validEnd = numberInSegment - discardNumber;

        if validEnd <= validStart
            error('区间过短，请减小discardRatio。');
        end

        segmentPower = segmentPower(validStart:validEnd);

        % 使用中位数比均值更抗异常脉冲
        meanPower(m,q) = median(segmentPower);
    end
end

%% =========================================================
%  9. 四相位差分恢复每个阵元的复数交叉项
%% =========================================================
P0   = meanPower(:,1);
P90  = meanPower(:,2);
P180 = meanPower(:,3);
P270 = meanPower(:,4);

realPart = (P0 - P180)/4;
imagPart = (P90 - P270)/4;

crossTerm = realPart + 1j*imagPart;

%% =========================================================
%  10. 本振通道校准
%
%  betaCalibration包含每个阵元本振链路的复增益。
%  初次测试可以全部设为1，但实际测角前最好标定。
%% =========================================================
betaCalibration = ones(4,1);

% 如果测得各通道幅相误差，可写成：
%
% betaCalibration = [
%     1.00*exp(1j*deg2rad(0));
%     0.96*exp(1j*deg2rad(8));
%     1.03*exp(1j*deg2rad(-5));
%     0.98*exp(1j*deg2rad(4))
% ];

estimatedArraySignal = crossTerm ./ conj(betaCalibration);

% 归一化，去掉未知的整体幅度
if norm(estimatedArraySignal) < eps
    error('恢复出的阵列信号接近零，请检查参考信号功率和相位切换。');
end

estimatedArraySignal = estimatedArraySignal ...
                       / norm(estimatedArraySignal);

%% =========================================================
%  11. 二维角度网格搜索
%
%  theta：离开阵列法向的夹角
%  phi：阵列平面内的方位角
%% =========================================================
thetaGrid = -60:0.5:60;
phiGrid   = -90:0.5:90;

spatialSpectrum = zeros(length(thetaGrid), ...
                        length(phiGrid));

for it = 1:length(thetaGrid)

    theta = deg2rad(thetaGrid(it));

    for ip = 1:length(phiGrid)

        phi = deg2rad(phiGrid(ip));

        ux = sin(theta)*cos(phi);
        uy = sin(theta)*sin(phi);

        phaseTerm = -k0 * ...
            (elementPosition(:,1)*ux + ...
             elementPosition(:,2)*uy);

        steeringVector = exp(1j*phaseTerm);
        steeringVector = steeringVector/norm(steeringVector);

        % 归一化匹配波束形成谱
        spatialSpectrum(it,ip) = ...
            abs(steeringVector' * estimatedArraySignal)^2;
    end
end

%% 找到最大值
[maxValue,maxLinearIndex] = max(spatialSpectrum(:));

[thetaIndex,phiIndex] = ind2sub( ...
    size(spatialSpectrum),maxLinearIndex);

estimatedTheta = thetaGrid(thetaIndex);
estimatedPhi = phiGrid(phiIndex);

fprintf('\n================ 测角结果 ================\n');
fprintf('估计俯仰参数 theta = %.2f°\n', estimatedTheta);
fprintf('估计方位角 phi     = %.2f°\n', estimatedPhi);
fprintf('最大归一化匹配值   = %.4f\n', maxValue);

%% =========================================================
%  12. 绘制二维角度谱
%% =========================================================
figure('Position',[180,100,950,650]);

imagesc(phiGrid,thetaGrid, ...
    10*log10(spatialSpectrum/max(spatialSpectrum(:)) + eps));

axis xy;
colorbar;

xlabel('方位角 \phi (°)');
ylabel('俯仰参数 \theta (°)');
title('四阵元相干能量测量二维角度谱');

hold on;

plot(estimatedPhi,estimatedTheta,'rx', ...
    'MarkerSize',12, ...
    'LineWidth',2);

hold off;

caxis([-20,0]);

%% =========================================================
%  13. 显示恢复出的阵元相对相位
%% =========================================================
relativePhase = angle(estimatedArraySignal ...
                      / estimatedArraySignal(1));

relativePhaseDegree = rad2deg(relativePhase);

fprintf('\n相对于左上阵元的相位：\n');

for m = 1:numberOfElements
    fprintf('%s：%.2f°\n', ...
        elementNames{m},relativePhaseDegree(m));
end

figure('Position',[300,180,800,450]);

stem(1:4,relativePhaseDegree, ...
    'LineWidth',1.5);

grid on;
box on;

xticks(1:4);
xticklabels(elementNames);

ylabel('相对于左上阵元的相位 (°)');
title('四阵元恢复相位');