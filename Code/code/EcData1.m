% 这个代码用于读取相应的矩阵
% 我们先给出计算的具体步骤
% a. 一个大循环读取所有的文件
% b. 针对于某一个具体的文件，我们首先读取它的第三列和第四列
% c. 读取到的第三列和第四列，把字符换成数字，并且把单位进行替换
% d. 按照新的排序规则，对其进行排序

clear all;
clc;


m=20;
G=zeros(m,m);

for m1=1:m
    file_path = "R110.21k_"+num2str(m1)+".csv";
    
    for m2=1:m
        fileID = fopen(file_path, 'r');
        
        % 读取第一行
        line_number = m2; % 例如，要读取第一行，将此值设为 1

        % 循环读取直到达到目标行
        for d = 1:line_number
            line = fgets(fileID);
            if ~ischar(line)
                error('到达文件结尾而无法读取指定行');
            end
        end
        
        % 关闭文件
        fclose(fileID);
        
        % 使用逗号分割第一行的内容
        C = strsplit(line, ',');
        
        % 获取第一个字段的内容
        first_field_str = strtrim(C{2}); % 移除首尾空格
        numeric_value_str = regexp(first_field_str, '[-+]?[0-9]*\.?[0-9]+', 'match');
        if ~isempty(numeric_value_str)
            numeric_value = str2double(numeric_value_str); % 转换为数字
            
            % 如果数字后面跟着 "u"，则乘以 0.001
            if endsWith(first_field_str, 'u')
                numeric_value = numeric_value * 0.000001;
            end
            
            if endsWith(first_field_str, 'm')
                numeric_value = numeric_value * 0.001;
            end
            
            if endsWith(first_field_str, 'n')
                numeric_value = numeric_value * (10^(-9));
            end
            
            if endsWith(first_field_str, 'k')
                numeric_value = numeric_value * (10^(3));
            end

            % 显示结果
            if size(numeric_value,2)==1
                G(m1,m2)=numeric_value;
            else
                G(m1,m2)=0;
            end
        end
        
        first_field_str = strtrim(C{3}); % 移除首尾空格
        numeric_value_str = regexp(first_field_str, '[-+]?[0-9]*\.?[0-9]+', 'match');
        if ~isempty(numeric_value_str)
            numeric_value = str2double(numeric_value_str); % 转换为数字
            G(m1,m2)=G(m1,m2)*exp(1i*numeric_value*2*pi/360);
%             1i*numeric_value
%             exp(1i*numeric_value*2*pi/360)
        end
        
    end
end

% 接下来按照G矩阵行列交换

% 原始数据
data =linspace(1,m,m);

% 使用冒泡排序按照自定义的比较规则排序
n = length(data);
for i = 1:n
    for j = 1:n-1
        if custom_compare(data(j+1), data(j))
            % 交换位置
            temp = data(j+1);
            data(j+1) = data(j);
            data(j) = temp;
        end
    end
end
[a,data]=sort(data);
G1=G(:,data);

% 频率
w=2*pi*5*10^5;
Cx=10*10^(-9);
Cy=5*10^(-9);
L1=100*10^(-6);


H=inv(G1); % 求逆
H=H-diag((1i*w*Cx-1/(1i*w*L1))*ones(1,m)); %减去对角元
H=H/(1i*w*Cx); %重新标定系数
H=H+diag(ones(1,m)); % 减去对角元1
H1=abs(H);
          
%对角化
[Ev,E]=eig(H,'vector');

subplot(1,2,1)
plot(real(E),imag(E),'o')
% axis([-5,5,-5,5])

subplot(1,2,2)
for m1=1:m
    plot(Ev(:,m1).*conj(Ev(:,m1))) 
    hold on;
end

figure()
surf(G1.*conj(G1))


%这个是比较函数
function result = custom_compare(x, y)
    str_x = num2str(x);
    str_y = num2str(y);
    len_x = length(str_x);
    len_y = length(str_y);
    
    % 按位比较
    for i = 1:min(len_x, len_y)
        if str_x(i) ~= str_y(i)
            result = str_x(i) < str_y(i);
            return;
        end
    end
    
    % 如果前面的位数都相同，则按照数字长度进行比较
    result = len_x < len_y;
end