clear all
% 读取文件内容
filename = 'Main.obj';
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
drag = [];

% 遍历每一行，查找包含Drag和Vol的行
for i = 1:length(lines)
    line = lines{i};
    
    % 检查是否包含Drag数据
    if contains(line, 'Drag = ') || contains(line, 'Vol = ')
        % 提取Drag值
        drag_match = regexp(line, 'Drag = ([\d\.]+)', 'tokens');
        if ~isempty(drag_match)
            drag_value = str2double(drag_match{1}{1});
            drag = [drag; drag_value];
        end
        
        % 提取Volume值
        vol_match = regexp(line, 'Vol = ([\d\.]+)', 'tokens');
        if ~isempty(vol_match)
            vol_value = str2double(vol_match{1}{1});
            Volume = [Volume; vol_value];
        end
    end
end

% 显示结果
disp('提取完成！数据统计:');
disp(['共提取 ', num2str(length(drag)), ' 个Drag值']);
disp(['共提取 ', num2str(length(Volume)), ' 个Volume值']);
Num = length(drag);
% 可选：保存到MAT文件
save('drag_volume_data.mat', 'drag', 'Volume');
disp('数据已保存到 drag_volume_data.mat');

