clear ; clc
tic

% 定义搜索范围和步长
% h_values = linspace(2.5, 10, 75);
h_values =6:10;
S_values = 1:100;

% 初始化最小值和对应的h、S
min_EC = inf;
best_h = NaN;
best_S = NaN;

% 初始化存储结果的矩阵
EC_matrix = zeros(length(h_values), length(S_values));

% 计算 EC 值的三维图并搜索最小值
for i = 1:length(h_values)
    for j = 1:length(S_values)
        h = h_values(i);
        S = S_values(j);
        z = [h, S];
        EC = Obj_funnew(z);
        EC_matrix(i, j) = EC;
        
        % 检查是否找到更小的 EC
        if EC < min_EC
            min_EC = EC;
            best_h = h;
            best_S = S;
        end
    end
end

% 创建网格数据
[S_grid, h_grid] = meshgrid(S_values, h_values);

% 绘制三维图
figure;
surf(S_grid, h_grid, EC_matrix);
title('EC vs h and S');
xlabel('S');
ylabel('h');
zlabel('EC');

% 显示最小 EC 及对应的 h 和 S
fprintf('最小的EC值：%f\n', min_EC);
fprintf('最小EC对应的h：%f\n', best_h);
fprintf('最小EC对应的S：%d\n', best_S);
toc