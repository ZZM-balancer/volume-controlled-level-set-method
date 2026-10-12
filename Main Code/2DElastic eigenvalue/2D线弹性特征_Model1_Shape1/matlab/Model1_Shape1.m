clear all
% 读取文件内容
filename = 'MainEigen.obj';
fid = fopen(filename, 'r');
if fid == -1
    error('无法打开文件: %s', filename);
end

% 读取所有行
lines = textscan(fid, '%s', 'Delimiter', '\n');
lines = lines{1};
fclose(fid);

% 初始化存储数组
Volume = [];
lambdaT = [];

% 遍历每一行，查找包含Drag和Vol的行
for i = 1:length(lines)
    line = lines{i};
    
    % 检查是否包含Drag数据
    if contains(line, 'lambda1 = ') || contains(line, 'Vol = ')
        % 提取Drag值
        lambda_match = regexp(line, 'lambda1 = ([\d\.]+)', 'tokens');
        if ~isempty(lambda_match)
            lambda_value = str2double(lambda_match{1}{1});
            lambdaT = [lambdaT; lambda_value];
        end
        
        % 提取Volume值
        vol_match = regexp(line, 'Vol = ([\d\.]+)', 'tokens');
        if ~isempty(vol_match)
            vol_value = str2double(vol_match{1}{1});
            Volume = [Volume; vol_value];
        end
    end
end
for i = 1:length(Volume)
    lambda(i) = lambdaT(4*i-3);
end
% 显示结果
disp('提取完成！数据统计:');
disp(['共提取 ', num2str(length(lambda)), ' 个Drag值']);
disp(['共提取 ', num2str(length(Volume)), ' 个Volume值']);
Num = length(Volume);
% 可选：保存到MAT文件
save('lambda_volume_data.mat', 'lambda', 'Volume');
disp('数据已保存到 lambda_volume_data.mat');

