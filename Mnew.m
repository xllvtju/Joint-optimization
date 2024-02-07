clear ; clc
tic

% ����������Χ�Ͳ���
% h_values = linspace(2.5, 10, 75);
h_values =6:10;
S_values = 1:100;

% ��ʼ����Сֵ�Ͷ�Ӧ��h��S
min_EC = inf;
best_h = NaN;
best_S = NaN;

% ��ʼ���洢����ľ���
EC_matrix = zeros(length(h_values), length(S_values));

% ���� EC ֵ����άͼ��������Сֵ
for i = 1:length(h_values)
    for j = 1:length(S_values)
        h = h_values(i);
        S = S_values(j);
        z = [h, S];
        EC = Obj_funnew(z);
        EC_matrix(i, j) = EC;
        
        % ����Ƿ��ҵ���С�� EC
        if EC < min_EC
            min_EC = EC;
            best_h = h;
            best_S = S;
        end
    end
end

% ������������
[S_grid, h_grid] = meshgrid(S_values, h_values);

% ������άͼ
figure;
surf(S_grid, h_grid, EC_matrix);
title('EC vs h and S');
xlabel('S');
ylabel('h');
zlabel('EC');

% ��ʾ��С EC ����Ӧ�� h �� S
fprintf('��С��ECֵ��%f\n', min_EC);
fprintf('��СEC��Ӧ��h��%f\n', best_h);
fprintf('��СEC��Ӧ��S��%d\n', best_S);
toc