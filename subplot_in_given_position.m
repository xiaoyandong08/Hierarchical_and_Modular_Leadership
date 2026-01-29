function sub_pos = subplot_in_given_position(outerPos,rows,cols)
%outerPos = [0.1, 0.05, 0.8, 0.5]; % figure 中的全局位置

% 分割为 2 行 3 列
% rows = 2;
% cols = 3;

% 计算每个子图的宽度和高度
dx = outerPos(3) / cols; % 每个子图的宽度
dy = outerPos(4) / rows; % 每个子图的高度


for r = 1:rows
    for c = 1:cols
        % 子图编号
        idx = (r - 1)*cols + c;
        
        % 计算当前子图的左下角坐标
        left = outerPos(1) + (c - 1)*dx;
        bottom = outerPos(2) + outerPos(4) - r*dy; % 注意从上往下排列
        
        % 创建 axes
        sub_pos(idx,:) = [left, bottom, dx, dy];
        
    end
end