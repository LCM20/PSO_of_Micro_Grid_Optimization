%% debug
% 设置导入选项并导入数据
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
kW = tbl.kW;
kW1 = tbl.kW1;
kW2 = tbl.kW2;
clear opts tbl
%使用前面的15分钟优化的程序单独加出来一天96个时段的总花费。应该与Optimization_Allday的结果相同。
cost = 0;
for time = 1:10
           N = 50;                         % 初始种群个数，所以应该有50个的以下三个变量。
            G = zeros(1,N); %电网输入的电能
            Gw = zeros(1,N);% 风力
            Gs = zeros(1,N);% 太阳能
            ger = 500;                      % 最大迭代次数  
            timeof96 = time;                  %运行表中第几个数据

            Gwlimit = zeros(1, 2);              %这里其实每个时间段，都有一个限制。需要回头改。  
            Gwlimit(1,2) = kW1(timeof96);             % 第一列是零，第二列是上限。
            Gslimit = zeros(1, 2);                
            Gslimit(1,2) = kW2(timeof96);

            vlimit = [-5, 5];               % 设置速度限制
            Gwv = vlimit(2)*rand(1,N);                  % 初始种群的速度
            Gsv = vlimit(2)*rand(1,N);                  % 初始种群的速度
            w = 0.8;                        % 惯性权重。%这里注意，看上面！我的速度设置的是一个0到1的随机数，但是实际上我的位置的变动的尺度可能是几百那么大，如此一来我的速度就非常小几乎可以忽略不记了。因此如果我想要达到一种比较好的效果，那么我就需要在这里设置一个系数。
            c1 = 0.1;                       % 自我学习因子
            c2 = 0.1;                       % 群体学习因子
            for i = 1:N
                Gw(1,i) = Gwlimit(1, 1) + (Gwlimit(1, 2) - Gwlimit(1, 1)) .* rand(1);%初始种群的位置
                Gs(1,i) = Gslimit(1, 1) + (Gslimit(1, 2) - Gslimit(1, 1)) .* rand(1);%初始种群的位置
            end

            Gwmx = Gw;                          % 每个个体的历史最佳位置
            Gsmx = Gs;                          % 每个个体的历史最佳位置
            Gwmxs = zeros(1);                % 种群的历史最佳位置，复数最后加s
            Gsmxs = zeros(1);                % 种群的历史最佳位置，复数最后加s

            %一个G、Gw、Gs对应一个个体，其三者可以计算出来一个适应度。
            G15minmy = inf(1,N);             % 每个个体的历史最佳适应度。
            G15minmys = inf(1);            % 种群历史最佳适应度

            % 群体更新
            iter = 1;
            record = zeros(1,ger);          % 记录器
            Gwrecord = repmat(Gw,ger,1);
            Gsrecord = repmat(Gs,ger,1);
            Grecord = repmat(G,ger,1);
            G15min = Gwrecord;              % G15min 行是次数，列是个体。

            GPrice = zeros(96,2);%第一列是售电电价，第二列是购电电价。
            GPrice(1:28,1) = 0.22;
            GPrice(1:28,2) = 0.25;
            GPrice(29:40,1) = 0.42;
            GPrice(29:40,2) = 0.53;
            GPrice(41:60,1) = 0.65;
            GPrice(41:60,2) = 0.82;
            GPrice(61:72,1) = 0.42;
            GPrice(61:72,2) = 0.53;
            GPrice(73:84,1) = 0.65;
            GPrice(73:84,2) = 0.82;
            GPrice(85:96,1) = 0.42;
            GPrice(85:96,2) = 0.53;


            while iter <= ger
                for i = 1:N
                    G(1,i) = kW(timeof96) - Gw(1,i) - Gs(1,i);
                end

                % G15min 行是次数，列是个体。
                for i = 1:N %G电网电量，为正为购入，为负则为卖出。
                    G15min(iter,i) = G(1,i).*(GPrice(timeof96,2).*(G(1,i)>0)+GPrice(timeof96,1).*(G(1,i)<=0)) + Gw(1,i)*0.52 + Gs(1,i)*0.75;
                end


                 for j = 1:N      
                    if G15minmy(1,j) > G15min(iter,j) % 这里是不是应该取最小值？？？
                        Gwmx(1,j) = Gw(1,j);
                        Gsmx(1,j) = Gs(1,j);   % 更新个体历史最佳位置
                        G15minmy(1,j) = G15min(iter,j);     % 更新个体历史最佳适应度
                    end 
                 end


                if G15minmys > min(G15min(iter,:))
                    [G15minmys,nmax] = min(G15min(iter,:)); % 更新群体历史最佳适应度
                    Gwmxs = Gw(1,nmax);      % 更新群体历史最佳位置
                    Gsmxs = Gs(1,nmax);      % 更新群体历史最佳位置
                end

                 Gwv = Gwv*w + c1*rand*(Gwmx - Gw) + c2*rand*(Gwmxs - Gw); 
                 Gsv = Gsv*w + c1*rand*(Gsmx - Gs) + c2*rand*(Gsmxs - Gs); 
                 %v = v * w + c1 * rand * (xm - x) + c2 * rand * (repmat(ym, N, 1) - x);% 速度更新
                 % 边界速度处理

                 Gwv(Gwv > vlimit(2)) = vlimit(2);
                 Gwv(Gwv < vlimit(1)) = vlimit(1);
                 Gsv(Gsv < vlimit(1)) = vlimit(1);
                 Gsv(Gsv > vlimit(2)) = vlimit(2);

                 Gwrecord(iter,:) = Gw; %行是个体的数值，列是第几次运算。
                 Gsrecord(iter,:) = Gs;
                 Grecord(iter,:) = G;
                 Gw = Gw + Gwv;% 位置更新
                 Gs = Gs + Gsv;% 位置更新
                 % 边界位置处理
                 Gw(Gw > Gwlimit(1,2)) =  Gwlimit(1,2);
                 Gs(Gs > Gslimit(1,2)) =  Gslimit(1,2);
                 Gw(Gw < Gwlimit(1,1)) =  Gwlimit(1,1);
                 Gs(Gs < Gslimit(1,1)) =  Gslimit(1,1);

                 record(1,iter) = G15minmys;

                 iter = iter+1;
            end 
            cost = cost + record(1,end);
end
cost/4