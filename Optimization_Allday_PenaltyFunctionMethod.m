%% 2020.12.4 考虑加入罚函数。详细安排都写在github的说明里了。
%% 这里用到的充电放电都是15分钟内的平均功率。
%% 设置导入选项并导入数据
clear;
data;
LoaPow = Pdata(:,2);
WinPowMax = Pdata(:,3);
SolPowMax = Pdata(:,4);
BatPowMax = 300*0.2/4;

% 写总的粒子群算法
MMMax = 99*2; % 这个是惩罚因子
h=waitbar(0,'Please wait');
N = 99;                  % 初始种群个数，所以应该有50个的以下三个变量。
GriPow96N = zeros(96,N); %电网输入的电能
ger = 999*3;                      % 最大迭代次数  

WinPowLimit962 = zeros(96, 2);              %这里其实每个时间段，都有一个限制。需要回头改。  
WinPowLimit962(:,2) = WinPowMax;             % 第一列是零，第二列是上限。
SolPowLimit962 = zeros(96, 2);                
SolPowLimit962(:,2) = SolPowMax;
BatPowLimit962 = zeros(96, 2);

%Parl = 56;
BatPowLimit962(:,2) = 300*0.2/4;          %蓄电池最大充放电功率
BatPowLimit962(:,1) = -300*0.2/4;          %蓄电池最大充放电功率

Vlimit12 = [-0.005, 0.005];               % 设置速度限制
XWinPowV96N = Vlimit12(2)*rand(96,N);                  % 初始种群的速度
XSolPowV96N = Vlimit12(2)*rand(96,N);                  % 初始种群的速度
XBatPowV96N = Vlimit12(2)*rand(96,N);              % 初始种群的速度

w = 0.8;                        % 惯性权重。%这里注意，看上面！我的速度设置的是一个0到1的随机数，但是实际上我的位置的变动的尺度可能是几百那么大，如此一来我的速度就非常小几乎可以忽略不记了。因此如果我想要达到一种比较好的效果，那么我就需要在这里设置一个系数。
c1 = 2;                       % 自我学习因子
c2 = 2;                       % 群体学习因子

XWinPow96N = rand(96,N);
XSolPow96N = rand(96,N);
XBatPow96N = rand(96,N);

WinPow96N = XWinPow96N .* WinPowMax;
SolPow96N = XSolPow96N .* SolPowMax;
BatPow96N = XBatPow96N .* BatPowMax * 2 - BatPowMax;



XWinPowMX96N = XWinPow96N;                   % 每个个体的历史最佳位置
XSolPowMX96N = XWinPow96N;                   % 每个个体的历史最佳位置
XBatPowMX96N = XWinPow96N;                   % 每个个体的历史最佳位置
XWinPowMXs961 = zeros(96,1);                % 种群的历史最佳位置，复数最后加s
XSolPowMXs961 = zeros(96,1);                % 种群的历史最佳位置，复数最后加s
XBatPowMXs961 = zeros(96,1);                % 种群的历史最佳位置，复数最后加s

%一个G、Gw、Gs对应一个个体，其三者可以计算出来一个适应度。
C15MY96N = inf(96,N);           % 每个个体的历史最佳适应度。
C15MYs961 = inf(96,1);          % 种群历史最佳适应度
CAlldayMY1N = inf(1,N);         % 以一整天的成本为适应度
CAlldayMYs = inf;
CAllday1N = CAlldayMY1N;          % Cost of all day.

% 群体更新
iter = 1;
RecordCost1ger = zeros(1,ger);          % 记录器
RecordCostAndPenalty1ger = zeros(1,ger);          % 记录器
Cost15_96N = zeros(96,N);              % 15分钟花费 行是次数，列是个体。
temp90 = zeros(96,N);
for i = 1   % 电网价格
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
end

figure;
while iter <= ger
    WinPow96N = XWinPow96N .* WinPowMax;
    SolPow96N = XSolPow96N .* SolPowMax;
    BatPow96N = XBatPow96N .* BatPowMax * 2 - BatPowMax;

    for i = 1:N
        GriPow96N(:,i) = LoaPow + BatPow96N(:,i) - WinPow96N(:,i) - SolPow96N(:,i);  % 
    end
    % G15min 行是次数，列是个体。
    for i = 1:N %G电网电量，为正为购入，为负则为卖出。
        temp125 = BatPow96N(:,i);
        Cost15_96N(:,i) = ( GriPow96N(:,i).*(GriPrice962(:,2).*(GriPow96N(:,i)>0)+GriPrice962(:,1).*(GriPow96N(:,i)<=0))...
            + WinPow96N(:,i)*0.52 + SolPow96N(:,i)*0.75)/4 - 0.2 * temp125 .*(temp125 < 0) ; %  注意这里应该是负号。
        
        % 下面这一行用来最终输出结果了，所以很重要。
        temp90(:,i) = ( GriPow96N(:,i).*(GriPrice962(:,2).*(GriPow96N(:,i)>0)+GriPrice962(:,1).*(GriPow96N(:,i)<=0))...
            + WinPow96N(:,i)*0.52 + SolPow96N(:,i)*0.75)/4  - 0.2 * temp125 .*(temp125 < 0) ... % 
            ;
    end

%     for i = 1:N %G电网电量，为正为购入，为负则为卖出。
%         Cost15_96Nger(:,i,iter) = ( GriPow96N(:,i).*(GriPrice962(:,2).*(GriPow96N(:,i)>0)+GriPrice962(:,1).*(GriPow96N(:,i)<=0))...
%             + WinPow96N(:,i)*0.52 + SolPow96N(:,i)*0.75)/4;
%     end

    CAllday1N = sum(Cost15_96N(41:end,:)) ...
            ...%      + MMMax*999999 * (max(0,  -( sum(BatPow96N(41:end,:)) + 300*0.6)               )).^2 ...
                  ...%+ MMMax*999999 * (max(0,   sum(BatPow96N(1:40,:))- 30                )).^2 ...
                ;

   for j = 1:N                                  %更新个体最小值
        if CAlldayMY1N(j) > CAllday1N(j)          % 这里是不是应该取最小值？？？
            XWinPowMX96N(:,j) = XWinPow96N(:,j);
            XSolPowMX96N(:,j) = XSolPow96N(:,j);   % 更新个体历史最佳位置
            XBatPowMX96N(:,j) = XBatPow96N(:,j);
            CAlldayMY1N(j) = CAllday1N(j);         % 更新个体历史最佳适应度     
        end 
   end 
    
    if CAlldayMYs > min(CAllday1N)
        [CAlldayMYs,nmax] = min(CAllday1N); % 更新群体历史最佳适应度
        XWinPowMXs961 = XWinPow96N(:,nmax);      % 更新群体历史最佳位置
        XSolPowMXs961 = XSolPow96N(:,nmax);      % 更新群体历史最佳位置
        XBatPowMXs961 = XBatPow96N(:,nmax);
    end

    XWinPowV96N = XWinPowV96N*w + c1*rand(96,N).*(XWinPowMX96N - XWinPow96N) + c2*rand(96,N).*(XWinPowMXs961 - XWinPow96N); 
    XSolPowV96N = XSolPowV96N*w + c1*rand(96,N).*(XSolPowMX96N - XSolPow96N) + c2*rand(96,N).*(XSolPowMXs961 - XSolPow96N); 
    XBatPowV96N = XBatPowV96N*w + c1*rand(96,N).*(XBatPowMX96N - XBatPow96N) + c2*rand(96,N).*(XBatPowMXs961 - XBatPow96N); 

     % 边界速度处理
     XWinPowV96N(XWinPowV96N < Vlimit12(1,1)) = Vlimit12(1,1);
     XWinPowV96N(XWinPowV96N > Vlimit12(1,2)) = Vlimit12(1,2);
     XSolPowV96N(XSolPowV96N < Vlimit12(1,1)) = Vlimit12(1,1);
     XSolPowV96N(XSolPowV96N > Vlimit12(1,2)) = Vlimit12(1,2);
     XBatPowV96N(XBatPowV96N < Vlimit12(1,1)) = Vlimit12(1,1);
     XBatPowV96N(XBatPowV96N > Vlimit12(1,2)) = Vlimit12(1,2);
     
     XWinPow96N = XWinPow96N + XWinPowV96N;% 位置更新
     XSolPow96N = XSolPow96N + XSolPowV96N;% 位置更新
     XBatPow96N = XBatPow96N + XBatPowV96N;% 位置更新
    
     XWinPow96N(XWinPow96N > 1) = 1;% 位置更新
     XWinPow96N(XWinPow96N < 0) = 0;% 位置更新
     XSolPow96N(XSolPow96N > 1) = 1;% 位置更新
     XSolPow96N(XSolPow96N < 0) = 0;% 位置更新
     XBatPow96N(XBatPow96N > 1) = 1;% 位置更新
     XBatPow96N(XBatPow96N < 0) = 0;% 位置更新

     RecordCost1ger(1,iter) = min(sum(temp90(41:end,:)));
     RecordCostAndPenalty1ger(1,iter) = min(CAllday1N);
     iter = iter+1;
     
     waitbar(iter/ger,h);
     
     % 画图
      if mod(iter,51) == 0

            surfx = subplot(2,3,1);
            surf(surfx,XWinPow96N(:,1:10),'Facecolor','interp','Edgecolor','interp')
            view(90,0)
            title(surfx,'X（风电功率）')
            
            surfx = subplot(2,3,2);
            surf(surfx,XSolPow96N(:,1:10),'Facecolor','interp','Edgecolor','interp')
            view(90,0)
            title(surfx,'X（光伏功率）')
            
            surfx = subplot(2,3,3);
            plot(surfx,linspace(0,24,96),XBatPow96N(:,1))
     
            
            title(surfx,'X（蓄电池功率）')
            
            pause(0.000001);
%         saveas(h,[num2str(iter),'.jpg'])
      end   
    
end
clear temp130;
clear temp125;
clear temp90;
delete(h);
figure;
subplot(2,1,1);
% 因为我们15分钟的时间，当成了一个算的，计算功率没有乘时间，所以这里除以4。
plot(RecordCost1ger); 
title("不加罚函数全天总花费");
subplot(2,1,2);
plot(RecordCostAndPenalty1ger);
title("有罚函数的全天总花费");
% subplot(2,2,2);
% plot(Grecord);
% title("电网电量");
% subplot(2,2,3);
% plot(Gwrecord);
% title("风力发电");
% subplot(2,2,4);
% plot(Gsrecord);
% title("太阳能发电");