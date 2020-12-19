%% 2020.12.4 考虑加入罚函数。详细安排都写在github的说明里了。
%% 这里用到的充电放电都是15分钟内的平均功率。
%% 设置导入选项并导入数据
clear;
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


% 写总的粒子群算法
temp21 = 0; % 这个变量记录的是，风电的功率在寻优过程中超出边界的次数。还挺多的，1038582这么多。因此加罚函数比较有必要。
MMMax = 99; % 这个是惩罚因子
h=waitbar(0,'Please wait');
N = 99;                  % 初始种群个数，所以应该有50个的以下三个变量。
GriPow96N = zeros(96,N); %电网输入的电能
WinPow96N = zeros(96,N);% 风力
SolPow96N = zeros(96,N);% 太阳能
BatPow96N = zeros(96,N);% 蓄电池电量。
ger = 99999;                      % 最大迭代次数  

WinPowLimit962 = zeros(96, 2);              %这里其实每个时间段，都有一个限制。需要回头改。  
WinPowLimit962(:,2) = WinPowMax;             % 第一列是零，第二列是上限。
SolPowLimit962 = zeros(96, 2);                
SolPowLimit962(:,2) = SolPowMax;
BatPowLimit962 = zeros(96, 2);
% 这里先假设前十个时间段是放电，然后后面十个个时间段是充电。


Parl = 40;
BatPowLimit962(1:Parl,2) = 300*0.2/4;          %蓄电池最大充放电功率
BatPowLimit962(1:Parl,1) = -300*0.2/4;          %蓄电池最大充放电功率

Vlimit12 = [-5, 5];               % 设置速度限制
WinPowV96N = Vlimit12(2)*rand(96,N);                  % 初始种群的速度
SolPowV96N = Vlimit12(2)*rand(96,N);                  % 初始种群的速度
BatPowV96N = Vlimit12(2)*rand(96,N) - Vlimit12(2)/2;              % 初始种群的速度

w = 0.8;                        % 惯性权重。%这里注意，看上面！我的速度设置的是一个0到1的随机数，但是实际上我的位置的变动的尺度可能是几百那么大，如此一来我的速度就非常小几乎可以忽略不记了。因此如果我想要达到一种比较好的效果，那么我就需要在这里设置一个系数。
c1 = 2;                       % 自我学习因子
c2 = 2;                       % 群体学习因子
% for i = 1:N
%     WinPow96N(:,i) = WinPowLimit962(:,1) + (WinPowLimit962(:,2) - WinPowLimit962(:,1)) .* rand(1);%初始种群的位置
%     SolPow96N(:,i) = SolPowLimit962(:,1) + (SolPowLimit962(:,2) - SolPowLimit962(:,1)) .* rand(1);%初始种群的位置
% end

% WinPow96N = 100*rand(96,N);%初始种群的位置
% SolPow96N = 100*rand(96,N);%初始种群的位置
% BatPow96N = 100*rand(96,N);%初始种群的位置

WinPow96N = zeros(96,N);%初始种群的位置
SolPow96N = zeros(96,N);%初始种群的位置
BatPow96N = zeros(96,N);%初始种群的位置
for i = 1:N
    
    BatPow96N(:,i) = BatPowLimit962(:,1) + (BatPowLimit962(:,2) - BatPowLimit962(:,1)) .* rand(1);%初始种群的位置
end

WinPowMX96N = WinPow96N;                   % 每个个体的历史最佳位置
SolPowMX96N = SolPow96N;                   % 每个个体的历史最佳位置
BatPowMX96N = BatPow96N;                   % 每个个体的历史最佳位置
WinPowMXs961 = zeros(96,1);                % 种群的历史最佳位置，复数最后加s
SolPowMXs961 = zeros(96,1);                % 种群的历史最佳位置，复数最后加s
BatPowMXs961 = zeros(96,1);                % 种群的历史最佳位置，复数最后加s

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
% WinPowRecord96Nger = zeros(96,N,ger);     % Gw = zeros(96,N);
% SolPowRecord96Nger = zeros(96,N,ger);
% GriPowRecord96Nger = zeros(96,N,ger);
% BatPowRecord96Nger = zeros(96,N,ger);
Cost15_96Nger = zeros(96,N,ger);              % 15分钟花费 行是次数，列是个体。
temp90 = zeros(96,N,ger);
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

% 只使用前几个时间段
while iter <= ger
    
%     if iter == 50
%        BatPow96N(1,1) = -3; 
%        BatPow96N(2,1) = -3; 
%     end
    for i = 1:N
        % 妙啊！BatPow有可能取值为负是电源，为正是负荷。发现不用标注充放电了！！！
        % 你仔细想啊，这里的蓄电池如果在某时间段小于零，说明在放电。这里与风电不一样，风电可以弃风弃光，不放出极限功率，
        % 而这里，如果小于零，那么意味着这个数值的电量是全都放完的。如果不想全都放完，我们就不在这里设置这么大的数值了。
        % 或者，其实蓄电池就是个功率可以为负的负载。
        GriPow96N(:,i) = LoaPow + BatPow96N(:,i) - WinPow96N(:,i) - SolPow96N(:,i);  % 
    end

    % G15min 行是次数，列是个体。
    for i = 1:N %G电网电量，为正为购入，为负则为卖出。
        temp125 = BatPow96N(:,i);
        Cost15_96Nger(:,i,iter) = ( GriPow96N(:,i).*(GriPrice962(:,2).*(GriPow96N(:,i)>0)+GriPrice962(:,1).*(GriPow96N(:,i)<=0))...
            + WinPow96N(:,i)*0.52 + SolPow96N(:,i)*0.75)/4 - 0.02 * temp125 .*(temp125 < 0)  ... %  注意这里应该是负号。
            + MMMax * (max(0,WinPow96N(:,i) - WinPowLimit962(:,2))).^2 ...
            + MMMax * (max(0,-WinPow96N(:,i) + WinPowLimit962(:,1))).^2 ...
            + MMMax * (max(0,SolPow96N(:,i) - SolPowLimit962(:,2))).^2 ...
            + MMMax * (max(0,-SolPow96N(:,i) + SolPowLimit962(:,1))).^2 ...
              + MMMax * (max(0,BatPow96N(:,i) - BatPowLimit962(:,2))).^2 ...
              + MMMax * (max(0,-BatPow96N(:,i) + BatPowLimit962(:,1))).^2 ...
        ;
        % 下面这一行用来最终输出结果了，所以很重要。
        temp90(:,i,iter) = ( GriPow96N(:,i).*(GriPrice962(:,2).*(GriPow96N(:,i)>0)+GriPrice962(:,1).*(GriPow96N(:,i)<=0))...
            + WinPow96N(:,i)*0.52 + SolPow96N(:,i)*0.75)/4  - 0.02 * temp125 .*(temp125 < 0) ... % 
            ;
    end
%     for i = 1:N %G电网电量，为正为购入，为负则为卖出。
%         Cost15_96Nger(:,i,iter) = ( GriPow96N(:,i).*(GriPrice962(:,2).*(GriPow96N(:,i)>0)+GriPrice962(:,1).*(GriPow96N(:,i)<=0))...
%             + WinPow96N(:,i)*0.52 + SolPow96N(:,i)*0.75)/4;
%     end

    CAllday1N = sum(Cost15_96Nger(1:Parl,:,iter)) ...
                  + MMMax*999999 * (max(0,  -( sum(BatPow96N(1:Parl,:)) )               )).^2 ...
                  + MMMax*999999 * (max(0,   sum(BatPow96N(1:28,:))-150                )).^2 ...
                ;

   for j = 1:N                                  %更新个体最小值
        if CAlldayMY1N(j) > CAllday1N(j)          % 这里是不是应该取最小值？？？
            WinPowMX96N(:,j) = WinPow96N(:,j);
            SolPowMX96N(:,j) = SolPow96N(:,j);   % 更新个体历史最佳位置
            BatPowMX96N(:,j) = BatPow96N(:,j);
            CAlldayMY1N(j) = CAllday1N(j);         % 更新个体历史最佳适应度     
        end 
   end 
      

    
    if CAlldayMYs > min(CAllday1N)
        [CAlldayMYs,nmax] = min(CAllday1N); % 更新群体历史最佳适应度
        WinPowMXs961 = WinPow96N(:,nmax);      % 更新群体历史最佳位置
        SolPowMXs961 = SolPow96N(:,nmax);      % 更新群体历史最佳位置
        BatPowMXs961 = BatPow96N(:,nmax);
    end

    

    WinPowV96N = WinPowV96N*w + c1*rand*(WinPowMX96N - WinPow96N) + c2*rand*(WinPowMXs961 - WinPow96N); 
    SolPowV96N = SolPowV96N*w + c1*rand*(SolPowMX96N - SolPow96N) + c2*rand*(SolPowMXs961 - SolPow96N); 
    BatPowV96N = BatPowV96N*w + c1*rand*(BatPowMX96N - BatPow96N) + c2*rand*(BatPowMXs961 - BatPow96N); 

     
     %v = v * w + c1 * rand * (xm - x) + c2 * rand * (repmat(ym, N, 1) - x);% 速度更新
     % 边界速度处理
     WinPowV96N(WinPowV96N < Vlimit12(1,1)) = Vlimit12(1,1);
     WinPowV96N(WinPowV96N > Vlimit12(1,2)) = Vlimit12(1,2);
     SolPowV96N(SolPowV96N < Vlimit12(1,1)) = Vlimit12(1,1);
     SolPowV96N(SolPowV96N > Vlimit12(1,2)) = Vlimit12(1,2);
     BatPowV96N(BatPowV96N < Vlimit12(1,1)) = Vlimit12(1,1);
     BatPowV96N(BatPowV96N > Vlimit12(1,2)) = Vlimit12(1,2);
     
     
%      WinPowRecord96Nger(:,:,iter) = WinPow96N; %行是个体的数值，列是第几次运算。
%      SolPowRecord96Nger(:,:,iter) = SolPow96N;
%      GriPowRecord96Nger(:,:,iter) = GriPow96N;
%      BatPowRecord96Nger(:,:,iter) = BatPow96N;
     
     
     WinPow96N = WinPow96N + WinPowV96N;% 位置更新
     SolPow96N = SolPow96N + SolPowV96N;% 位置更新
     BatPow96N = BatPow96N + BatPowV96N;% 位置更新
     
     for i = 1 : 96
        for j = 1 : N
            if (WinPow96N(i,j) > WinPowLimit962(i,2)) || (WinPow96N(i,j) < WinPowLimit962(i,1))
                temp21 = temp21 + 1;
            end
        end
     end % 用来发现到底是谁超过了边界，发现还挺多。
     
     % 边界位置处理，有了惩罚项，不用约束边界了。
%      temp130 = repmat(WinPowLimit962(:,2),1,N);
%      WinPow96N(WinPow96N > temp130) =  temp130(WinPow96N > temp130);
%      temp130 = repmat(SolPowLimit962(:,2),1,N);
%      SolPow96N(SolPow96N > temp130) =  temp130(SolPow96N > temp130);
%      temp130 = repmat(WinPowLimit962(:,1),1,N);
%      WinPow96N(WinPow96N < temp130) =  temp130(WinPow96N < temp130);
%      temp130 = repmat(SolPowLimit962(:,1),1,N);
%      SolPow96N(SolPow96N < temp130) =  temp130(SolPow96N < temp130);

     RecordCost1ger(1,iter) = min(sum(temp90(1:Parl,:,iter)));
     RecordCostAndPenalty1ger(1,iter) = min(CAllday1N);
     iter = iter+1;
     
     waitbar(iter/ger,h);
    
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
%% 连续运行多次
for temp203 = 1:20
    temp21 = 0; % 这个变量记录的是，风电的功率在寻优过程中超出边界的次数。还挺多的，1038582这么多。因此加罚函数比较有必要。
    MMMax = 99; % 这个是惩罚因子
    h=waitbar(0,'Please wait');
    N = 99;                  % 初始种群个数，所以应该有50个的以下三个变量。
    ger = 999;                      % 最大迭代次数  

    % 这里先假设前十个时间段是放电，然后后面十个个时间段是充电。

    Parl = 96;
    
    %一个G、Gw、Gs对应一个个体，其三者可以计算出来一个适应度。
    
    % 群体更新
    iter = 1;
    
    % 此处规定充放电的时间。直接按照序号来。一共96；
    
    % 只使用前几个时间段
    while iter <= ger

        for i = 1:N
            % 妙啊！BatPow有可能取值为负是电源，为正是负荷。发现不用标注充放电了！！！
            % 你仔细想啊，这里的蓄电池如果在某时间段小于零，说明在放电。这里与风电不一样，风电可以弃风弃光，不放出极限功率，
            % 而这里，如果小于零，那么意味着这个数值的电量是全都放完的。如果不想全都放完，我们就不在这里设置这么大的数值了。
            % 或者，其实蓄电池就是个功率可以为负的负载。
            GriPow96N(:,i) = LoaPow + BatPow96N(:,i) - WinPow96N(:,i) - SolPow96N(:,i);  % 
        end

        % G15min 行是次数，列是个体。
        for i = 1:N %G电网电量，为正为购入，为负则为卖出。
            temp125 = BatPow96N(:,i);
            Cost15_96Nger(:,i,iter) = ( GriPow96N(:,i).*(GriPrice962(:,2).*(GriPow96N(:,i)>0)+GriPrice962(:,1).*(GriPow96N(:,i)<=0))...
                + WinPow96N(:,i)*0.52 + SolPow96N(:,i)*0.75)/4 - 0.02 * temp125 .*(temp125 < 0) ... %  注意这里应该是负号。
                + MMMax * (max(0,WinPow96N(:,i) - WinPowLimit962(:,2))).^2 ...
                + MMMax * (max(0,-WinPow96N(:,i) + WinPowLimit962(:,1))).^2 ...
                + MMMax * (max(0,SolPow96N(:,i) - SolPowLimit962(:,2))).^2 ...
                + MMMax * (max(0,-SolPow96N(:,i) + SolPowLimit962(:,1))).^2 ...
              + MMMax * (max(0,BatPow96N(:,i) - BatPowLimit962(:,2))).^2 ...
              + MMMax * (max(0,-BatPow96N(:,i) + BatPowLimit962(:,1))).^2 ...
            ;
            % 下面这一行用来最终输出结果了，所以很重要。
            temp90(:,i,iter) = ( GriPow96N(:,i).*(GriPrice962(:,2).*(GriPow96N(:,i)>0)+GriPrice962(:,1).*(GriPow96N(:,i)<=0))...
                + WinPow96N(:,i)*0.52 + SolPow96N(:,i)*0.75)/4 - 0.02 * temp125 .*(temp125 < 0) ... % 
                ;
        end
    %     for i = 1:N %G电网电量，为正为购入，为负则为卖出。
    %         Cost15_96Nger(:,i,iter) = ( GriPow96N(:,i).*(GriPrice962(:,2).*(GriPow96N(:,i)>0)+GriPrice962(:,1).*(GriPow96N(:,i)<=0))...
    %             + WinPow96N(:,i)*0.52 + SolPow96N(:,i)*0.75)/4;
    %     end

        CAllday1N = sum(Cost15_96Nger(1:Parl,:,iter)) ...
                  + MMMax*9999 * (max(0,  -( sum(BatPow96N(1:Parl,:)) )               )).^2 ...
                    ;

       for j = 1:N                                  %更新个体最小值
            if CAlldayMY1N(j) > CAllday1N(j)          % 这里是不是应该取最小值？？？
                WinPowMX96N(:,j) = WinPow96N(:,j);
                SolPowMX96N(:,j) = SolPow96N(:,j);   % 更新个体历史最佳位置
                BatPowMX96N(:,j) = BatPow96N(:,j);
                CAlldayMY1N(j) = CAllday1N(j);         % 更新个体历史最佳适应度     
            end 
       end 



        if CAlldayMYs > min(CAllday1N)
            [CAlldayMYs,nmax] = min(CAllday1N); % 更新群体历史最佳适应度
            WinPowMXs961 = WinPow96N(:,nmax);      % 更新群体历史最佳位置
            SolPowMXs961 = SolPow96N(:,nmax);      % 更新群体历史最佳位置
            BatPowMXs961 = BatPow96N(:,nmax);
        end



        WinPowV96N = WinPowV96N*w + c1*rand*(WinPowMX96N - WinPow96N) + c2*rand*(WinPowMXs961 - WinPow96N) + rand; 
        SolPowV96N = SolPowV96N*w + c1*rand*(SolPowMX96N - SolPow96N) + c2*rand*(SolPowMXs961 - SolPow96N) + rand; 
        BatPowV96N = BatPowV96N*w + c1*rand*(BatPowMX96N - BatPow96N) + c2*rand*(BatPowMXs961 - BatPow96N) + rand; 


         %v = v * w + c1 * rand * (xm - x) + c2 * rand * (repmat(ym, N, 1) - x);% 速度更新
         % 边界速度处理
         WinPowV96N(WinPowV96N < Vlimit12(1,1)) = Vlimit12(1,1);
         WinPowV96N(WinPowV96N > Vlimit12(1,2)) = Vlimit12(1,2);
         SolPowV96N(SolPowV96N < Vlimit12(1,1)) = Vlimit12(1,1);
         SolPowV96N(SolPowV96N > Vlimit12(1,2)) = Vlimit12(1,2);
         BatPowV96N(BatPowV96N < Vlimit12(1,1)) = Vlimit12(1,1);
         BatPowV96N(BatPowV96N > Vlimit12(1,2)) = Vlimit12(1,2);


%          WinPowRecord96Nger(:,:,iter) = WinPow96N; %行是个体的数值，列是第几次运算。
%          SolPowRecord96Nger(:,:,iter) = SolPow96N;
%          GriPowRecord96Nger(:,:,iter) = GriPow96N;
%          BatPowRecord96Nger(:,:,iter) = BatPow96N;


         WinPow96N = WinPow96N + WinPowV96N;% 位置更新
         SolPow96N = SolPow96N + SolPowV96N;% 位置更新
         BatPow96N = BatPow96N + BatPowV96N;% 位置更新

         for i = 1 : 96
            for j = 1 : N
                if (WinPow96N(i,j) > WinPowLimit962(i,2)) || (WinPow96N(i,j) < WinPowLimit962(i,1))
                    temp21 = temp21 + 1;
                end
            end
         end % 用来发现到底是谁超过了边界，发现还挺多。

         % 边界位置处理，有了惩罚项，不用约束边界了。
    %      temp130 = repmat(WinPowLimit962(:,2),1,N);
    %      WinPow96N(WinPow96N > temp130) =  temp130(WinPow96N > temp130);
    %      temp130 = repmat(SolPowLimit962(:,2),1,N);
    %      SolPow96N(SolPow96N > temp130) =  temp130(SolPow96N > temp130);
    %      temp130 = repmat(WinPowLimit962(:,1),1,N);
    %      WinPow96N(WinPow96N < temp130) =  temp130(WinPow96N < temp130);
    %      temp130 = repmat(SolPowLimit962(:,1),1,N);
    %      SolPow96N(SolPow96N < temp130) =  temp130(SolPow96N < temp130);

         RecordCost1ger(1,iter) = min(sum(temp90(1:Parl,:,iter)));
         RecordCostAndPenalty1ger(1,iter) = min(CAllday1N);
         iter = iter+1;

         waitbar(iter/ger,h);

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
 
end