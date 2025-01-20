% AL region 
% generalized boundary condition
% eigenvalue change


L = 100;
omega = vpa((sqrt(5)-1)/2, 50);  % 使用 vpa 提高精度到 50 位
lambda = vpa(1.5, 50);  % 使用 vpa 提高精度到 50 位

% 高精度矩阵构造
H = diag(ones(1,L-1),-1) + diag(2*lambda*cos(2*pi*omega*linspace(0,L-1,L)));

% 设置 H(1, L) 元素为 10^(-4)
H(1,L) = vpa(0, 50);  % 使用 vpa 确保精度

% H(1,L) = vpa(0,50); 
% 为了确保矩阵 H 是精确的，使用 vpa 函数转换成高精度类型
H = vpa(H, 50);  

% 计算特征值和特征向量
[Ev, D] = eig(H);  % 不使用 'vector'，直接计算特征值和特征向量


E1=diag(double(D));



% 设置 m_all 的值
m_all = [50, 70, 90];

% 初始化数据矩阵
data = zeros(L, 3);

Em=[2.4457,-2.9009,1.2715];

n = 1;
for m = m_all
    omgea=(sqrt(5)-1)/2;
    [~, idx] = min(abs(E1 -Em(n) ));
    % 计算归一化的特征向量的平方模
    
    ta = sum(Ev(:,idx).*conj(Ev(:,idx)));
    
    % 使用 semilogy 绘制特征向量的平方模
    semilogy(Ev(:,idx).*conj(Ev(:,idx))/ta, 'linewidth', 3)
    hold on;
    
    % 存储数据
    data(:,n) = Ev(:,idx).*conj(Ev(:,idx));
    n = n + 1;
end

k= 2*log(lambda);

X=linspace(10,L-10,L-10);
Y=exp(-k*X);
Y=Y/sqrt(sum(Y.*Y));

Y=Y/(10^10);
semilogy(X,Y,'k--','linewidth',3)

xlim([0,L])
ylim([10^(-50),1])

xlabel('$j$','interpreter','latex')
ylabel('$|\psi|^2$','interpreter','latex')

% 将 data 转换为数值类型
numericData = double(data);

% 保存 data 到 CSV 文件
save('psi2_AL_OBC.mat', 'numericData');
