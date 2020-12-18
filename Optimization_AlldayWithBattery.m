%% 2020.12.4 考虑加入罚函数。详细安排都写在github的说明里了。
%% 这里用到的充电放电都是15分钟内的平均功率。
%% 设置导入选项并导入数据
for i = 1 % 这里的for不循环，只是为了折叠代码。为什么不用if，是因为if没有折叠功能！？
    clear;
    opts = spreadsheetImportOptions("NumVariables", 3);
    % 指定工作表和范围
    opts.Sheet = "Sheet1";
    opts.DataRange = "B2:D97";
    % 指定列名称和类型
    opts.VariableNames = ["kW", "kW1", "kW2"];
    opts.VariableTypes = ["double", "double", "double"];
    % 导入数据
    tbl = readtable("C:\My_Files\课程\电力系统优化方法\优化作业\A题数据.xlsx", opts, "UseExcel", false);
    LoaPow = tbl.kW;
    WinPowMax = tbl.kW1;
    SolPowMax = tbl.kW2;
    clear opts tbl
end
%% 写总的粒子群算法
temp21 = 0; % 这个变量记录的是，风电的功率在寻优过程中超出边界的次数。还挺多的，1038582这么多。因此加罚函数比较有必要。
MMMax = 999; % 这个是惩罚因子
h=waitbar(0,'Please wait');

N = 50;                         % 初始种群个数，所以应该有50个的以下三个变量。
GriPow96N = zeros(96,N); %电网输入的电能
WinPow96N = zeros(96,N);% 风力
SolPow96N = zeros(96,N);% 太阳能
BatPow96N = zeros(96,N);% 蓄电池电量。
ger = 500;                      % 最大迭代次数  

WinPowLimit962 = zeros(96, 2);              %这里其实每个时间段，都有一个限制。需要回头改。  
WinPowLimit962(:,2) = WinPowMax;             % 第一列是零，第二列是上限。
SolPowLimit962 = zeros(96, 2);                
SolPowLimit962(:,2) = SolPowMax;
BatPowLimit962 = 300*0.2;          %蓄电池最大充放电功率

Vlimit12 = [-5, 5];               % 设置速度限制
WinPowV96N = Vlimit12(2)*rand(96,N);                  % 初始种群的速度
SolPowV96N = Vlimit12(2)*rand(96,N);                  % 初始种群的速度
BatPowV96N = Vlimit12(2)*rand(96,N);              % 初始种群的速度

w = 0.8;                        % 惯性权重。%这里注意，看上面！我的速度设置的是一个0到1的随机数，但是实际上我的位置的变动的尺度可能是几百那么大，如此一来我的速度就非常小几乎可以忽略不记了。因此如果我想要达到一种比较好的效果，那么我就需要在这里设置一个系数。
c1 = 0.1;                       % 自我学习因子
c2 = 0.1;                       % 群体学习因子
% for i = 1:N
%     WinPow96N(:,i) = WinPowLimit962(:, 1) + (WinPowLimit962(:, 2) - WinPowLimit962(:, 1)) .* rand(1);%初始种群的位置
%     SolPow96N(:,i) = SolPowLimit962(:, 1) + (SolPowLimit962(:, 2) - SolPowLimit962(:, 1)) .* rand(1);%初始种群的位置
% end

WinPow96N = 100*rand(96,N);%初始种群的位置
SolPow96N = 100*rand(96,N);%初始种群的位置
BatPow96N = 100*rand(96,N);%初始种群的位置

WinPowMX96N = WinPow96N;                          % 每个个体的历史最佳位置
SolPowMX96N = SolPow96N;                          % 每个个体的历史最佳位置
BatPowMX96N = BatPow96N;
WinPowMXs961 = zeros(96,1);                % 种群的历史最佳位置，复数最后加s
SolPowMXs961 = zeros(96,1);                % 种群的历史最佳位置，复数最后加s
BatPowMXs961 = zeros(96,1);                % 种群的历史最佳位置，复数最后加s

%一个G、Gw、Gs对应一个个体，其三者可以计算出来一个适应度。
C15MY96N = inf(96,N);           % 每个个体的历史最佳适应度。
C15MYs961 = inf(96,1);            % 种群历史最佳适应度


% 群体更新
iter = 1;
Record1ger = zeros(1,ger);          % 记录器
WinPowRecord96Nger = zeros(96,N,ger);     % Gw = zeros(96,N);
SolPowRecord96Nger = zeros(96,N,ger);
GriPowRecord96Nger = zeros(96,N,ger);
BatPowRecord96Nger = zeros(96,N,ger);
Cost15_96Nger = zeros(96,N,ger);              % 15分钟花费 行是次数，列是个体。

GriPrice962 = zeros(96,2);%第一列是售电电价，第二列是购电电价。
GriPrice962(1:28,1) = 0.22;
GriPrice962(1:28,2) = 0.25;
GriPrice962(29:40,1) = 0.42;
GriPrice962(29:40,2) = 0.53;
GriPrice962(41:60,1) = 0.65;
GriPrice962(41:60,2) = 0.82;
GriPrice962(61:72,1) = 0.42;
GriPrice962(61:72,2) = 0.53;
GriPrice962(73:84,1) = 0.65;
GriPrice962(73:84,2) = 0.82;
GriPrice962(85:96,1) = 0.42;
GriPrice962(85:96,2) = 0.53;

% 此处规定充放电的时间。直接按照序号来。一共96；
DisAndCha163 = zeros(16,2); %第一列是开始时间，第二列是结束时间，第三列用来标记充电还是放电。充电为1放电为0
temp84 = 1;
temp85 = 6;
for i = 1:16
    DisAndCha163(i,1) = temp84;
    DisAndCha163(i,2) = temp85;
    temp84 = temp84 + 6;
    temp85 = temp85 + 6;
end
while iter <= ger
    for i = 1:N
        % 妙啊！BatPow有可能取值为负是电源，为正是负荷。发现不用标注充放电了！！！
        % 你仔细想啊，这里的蓄电池如果在某时间段小于零，说明在放电。这里与风电不一样，风电可以弃风弃光，不放出极限功率，
        % 而这里，如果小于零，那么意味着这个数值的电量是全都放完的。如果不想全都放完，我们就不在这里设置这么大的数值了。
        % 或者，其实蓄电池就是个功率可以为负的负载。
        GriPow96N(:,i) = LoaPow + BatPow96N(:,i) - WinPow96N(:,i) - SolPow96N(:,i);  
    end

    % G15min 行是次数，列是个体。
    for i = 1:N %G电网电量，为正为购入，为负则为卖出。
        Cost15_96Nger(:,i,iter) = ( GriPow96N(:,i).*(GriPrice962(:,2).*(GriPow96N(:,i)>0)+GriPrice962(:,1).*(GriPow96N(:,i)<=0))...
            + WinPow96N(:,i)*0.52 + SolPow96N(:,i)*0.75)/4 ...
            + MMMax * (max(0,WinPow96N(:,i) - WinPowLimit962(:,2))).^2 ...
            + MMMax * (max(0,-WinPow96N(:,i) + WinPowLimit962(:,1))).^2 ...
            + MMMax * (max(0,SolPow96N(:,i) - SolPowLimit962(:,2))).^2 ...
            + MMMax * (max(0,-SolPow96N(:,i) + SolPowLimit962(:,1))).^2 ;
        for j = 1:16
            Cost15_96Nger(:,i,iter);
        end
        
    end
    
    for i = 1:96
       for j = 1:N      
            if C15MY96N(i,j) > Cost15_96Nger(i,j,iter) % 这里是不是应该取最小值？？？
                WinPowMX96N(i,j) = WinPow96N(i,j);
                SolPowMX96N(i,j) = SolPow96N(i,j);   % 更新个体历史最佳位置
                C15MY96N(i,j) = Cost15_96Nger(i,j,iter);     % 更新个体历史最佳适应度
            end 
       end 
    end

    for i = 1:96
        if C15MYs961(i) > min(Cost15_96Nger(i,:,iter))
            [C15MYs961(i),nmax] = min(Cost15_96Nger(i,:,iter)); % 更新群体历史最佳适应度
            WinPowMXs961(i) = WinPow96N(i,nmax);      % 更新群体历史最佳位置
            SolPowMXs961(i) = SolPow96N(i,nmax);      % 更新群体历史最佳位置
        end
    end
    
    for i = 1:96
        WinPowV96N(i,:) = WinPowV96N(i,:)*w + c1*rand*(WinPowMX96N(i,:) - WinPow96N(i,:)) + c2*rand*(WinPowMXs961(i) - WinPow96N(i,:)); 
        SolPowV96N(i,:) = SolPowV96N(i,:)*w + c1*rand*(SolPowMX96N(i,:) - SolPow96N(i,:)) + c2*rand*(SolPowMXs961(i) - SolPow96N(i,:)); 
        BatPowV96N(i,:) = BatPowV96N(i,:)*w + c1*rand*(BatPowMX96N(i,:) - BatPow96N(i,:)) + c2*rand*(BatPowMXs961(i) - BatPow96N(i,:)); 
    end
     
     %v = v * w + c1 * rand * (xm - x) + c2 * rand * (repmat(ym, N, 1) - x);% 速度更新
     % 边界速度处理
     
     % 这里应该设置下限的，可以后面再加。因为上一个程序没加。
     WinPowV96N(WinPowV96N < Vlimit12(1,1)) = Vlimit12(1,1);
     WinPowV96N(WinPowV96N > Vlimit12(1,2)) = Vlimit12(1,2);
     SolPowV96N(SolPowV96N < Vlimit12(1,1)) = Vlimit12(1,1);
     SolPowV96N(SolPowV96N > Vlimit12(1,2)) = Vlimit12(1,2);
     
     WinPowRecord96Nger(:,:,iter) = WinPow96N; %行是个体的数值，列是第几次运算。
     SolPowRecord96Nger(:,:,iter) = SolPow96N;
     GriPowRecord96Nger(:,:,iter) = GriPow96N;
     WinPow96N = WinPow96N + WinPowV96N;% 位置更新
     SolPow96N = SolPow96N + SolPowV96N;% 位置更新
     % 边界位置处理
     temp130 = repmat(WinPowLimit962(:,2),1,N);
     WinPow96N(WinPow96N > temp130) =  temp130(WinPow96N > temp130);
     temp130 = repmat(SolPowLimit962(:,2),1,N);
     SolPow96N(SolPow96N > temp130) =  temp130(SolPow96N > temp130);
     temp130 = repmat(WinPowLimit962(:,1),1,N);
     WinPow96N(WinPow96N < temp130) =  temp130(WinPow96N < temp130);
     temp130 = repmat(SolPowLimit962(:,1),1,N);
     SolPow96N(SolPow96N < temp130) =  temp130(SolPow96N < temp130);

     Record1ger(1,iter) = sum(C15MYs961);

     iter = iter+1;
     
     waitbar(iter/ger,h);
end
clear temp130;
clear temp111;
delete(h);
figure;
subplot(2,2,1);
% 因为我们15分钟的时间，当成了一个算的，计算功率没有乘时间，所以这里除以4。
plot(Record1ger); 
title("全天总花费");
% subplot(2,2,2);
% plot(Grecord);
% title("电网电量");
% subplot(2,2,3);
% plot(Gwrecord);
% title("风力发电");
% subplot(2,2,4);
% plot(Gsrecord);
% title("太阳能发电");