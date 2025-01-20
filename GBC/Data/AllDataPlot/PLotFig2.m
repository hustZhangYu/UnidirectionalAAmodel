% This code, we plot the figure with all of the datas
clear all;
clc;

main_plot()



function main_plot()
    % 主程序设置图形窗口
    figure('Position', [0, 0, 1800, 1000]);
    
    % 设置子图布局：2行4列，共8个子图
    subplot(2, 4, 1);
    plot_subgraph_1();
    
    subplot(2, 4, 2);
    plot_subgraph_2();
    
    subplot(2, 4, 3);
    plot_subgraph_3();
    
    subplot(2, 4, 4);
    plot_subgraph_4();
    
    subplot(2, 4, 5);
    plot_subgraph_5();
    
    subplot(2, 4, 6);
    plot_subgraph_6();
    
    subplot(2, 4, 7);
    plot_subgraph_7();
    
    subplot(2, 4, 8);
    plot_subgraph_8();
    
    figure()
    plot_inset1()
    
    figure()
    plot_inset2()
    

    
end

% 子图1的绘制
function plot_subgraph_1()
    lambda=0.5;
    L=100;
    a=load('psi2_SE_OBC.mat');
    data=a.data;
    semilogy(data,'linewidth',3)

    hold on;
    k= 2*log(lambda);

    X=linspace(60,L-10,L-10);
    Y=exp(-k*X);
    Y=Y/sqrt(sum(Y.*Y));

    Y=Y/(10^10);
    semilogy(X,Y,'k--','linewidth',3)

    xlabel('$j$','interpreter','latex')
    ylabel('$|\psi_m(j)|^2$','interpreter','latex')
    lgd=legend('$m=50$','$m=70$','$m=90$','Predicted','interpreter','latex','location','northwest')
    set(lgd, 'Color', [1, 1, 1, 0.5]);
    xlim([0,L])
    ylim([10^(-50),1])
    

end
% 子图2的绘制
function plot_subgraph_2()
    lambda = 0.5;
    L = 100;
    a = load('psi2_SE_GBC.mat');
    data = a.data;
    semilogy(data, 'linewidth', 3)
    hold on;
    
    xlabel('$j$', 'interpreter', 'latex')
    ylabel('$|\psi_m(j)|^2$', 'interpreter', 'latex')
    lgd = legend('$m=50$', '$m=70$', '$m=90$', 'interpreter', 'latex', 'location', 'southeast');
    set(lgd, 'Color', [1, 1, 1, 0.5]);
    xlim([0, L])
    ylim([10^(-20), 1])
    

end
% 子图3的绘制
function plot_subgraph_3()

L = 100;
gamma_all = [0.0001, 1, 10000];
dataName = {'SFL_Spectrum_gamma_0.0001.mat','AL_Spectrum_gamma_0.0001.mat','SFL_Spectrum_gamma_1.0000.mat','AL_Spectrum_gamma_1.0000.mat','SFL_Spectrum_gamma_10000.0000.mat','AL_Spectrum_gamma_10000.0000.mat'};
% 定义不同的标记符号
markers = ['o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*']; 
Lines = {'-', '--', '-.'};
gamma2=[-4,0,4];
n = 1;
% 创建一个空的 legend 条目列表
legendEntries = {};

for k = 1:6
    marker = markers(n);
    a = load(dataName{n});
    E = a.E;
    if mod(k, 2) == 1
        plot(imag(E), real(E), marker, 'markersize', 3)
        legendEntries{end+1} = ['$\lambda=0.5,\eta =10^{ ', num2str(gamma2((k+1)/2)),'}$'];
        hold on;
    end
    
    % 为每一组数据添加图例项
    n = n + 1;
end

n=1;
for k = 1:6
    marker = markers(n);
    a = load(dataName{n});
    E = a.E;
    if mod(k, 2) == 0
        plot(imag(E), real(E), marker, 'markersize', 3)
        legendEntries{end+1} = ['$\lambda=1.5,\eta =10^{ ', num2str(gamma2(k/2)),'}$'];
    end
    
    % 为每一组数据添加图例项
    n = n + 1; 
end
view(270, 90)

lambda = 0.5;
a = -3.5:0.1:3.5;
b = -3:0.1:3;
data = zeros(length(a), length(b));

for z = 1:length(gamma_all)
    gamma1 = gamma_all(z);
    for m = 1:length(a)
        for n = 1:length(b)
            a1 = a(m);
            b1 = b(n);
            c = a1 + 1i*b1;
            r1 = a1^2 / (1 + (2*lambda)^2 / (4*exp(2*log(gamma1)/L)))^2 + b1^2 / (1 - (2*lambda)^2 / (4*exp(2*log(gamma1)/L)))^2;
            data(m, n) = r1;
        end
    end
    line = Lines{z};
    contour(b, a, real(data), [exp(2*log(gamma1)/L), exp(2*log(gamma1)/L)], line, 'linewidth', 1, 'color', 'k')
    hold on;
    
    % 为每一条contour曲线添加图例项
%     legendEntries{end+1} = ['$\eta=10^{ ', num2str(gamma2(z)),'}$'];
end

% 添加所有的legend条目
legend(legendEntries, 'Location', 'Best' ,'NumColumns', 2,'interpreter','latex');

xlim([-3,2])
xlabel('$\Im (E)$','interpreter','latex')
ylabel('$\Re (E)$','interpreter','latex')
end
% 子图4的绘制
function plot_subgraph_4()

% 加载第一个数据文件
a = load('Eigenstates_SFL_bd_gamma1e-04.mat');
Ev = a.Ev;
L = a.L;
color = [0, 0.4470, 0.7410];
for i = 1:L
    if i == 1  % 只在第一次绘制时添加图例
        z1=semilogy(linspace(1, L, L), Ev(:, i).*conj(Ev(:, i)), 'linewidth',3,'color', color);
    else
        semilogy(linspace(1, L, L), Ev(:, i).*conj(Ev(:, i)), 'color', color);
    end
    hold on;
end

% 加载第二个数据文件
a = load('Eigenstates_SFL_bd_gamma1e+00.mat');
Ev = a.Ev;
L = a.L;
color = [0.8500, 0.3250, 0.0980];
for i = 1:L
    if i == 1  % 只在第一次绘制时添加图例
        z2=semilogy(linspace(1, L, L), Ev(:, i).*conj(Ev(:, i)), 'linewidth',3, 'color', color);
    else
        semilogy(linspace(1, L, L), Ev(:, i).*conj(Ev(:, i)), 'color', color);
    end
    hold on;
end

% 加载第三个数据文件
a = load('Eigenstates_SFL_bd_gamma1e+04.mat');
Ev = a.Ev;
L = a.L;
color = [0.4660, 0.6740, 0.1880];
for i = 1:L
    if i == 1  % 只在第一次绘制时添加图例
        z3=semilogy(linspace(1, L, L), Ev(:, i).*conj(Ev(:, i)), 'linewidth',3, 'color', color);
    else
        semilogy(linspace(1, L, L), Ev(:, i).*conj(Ev(:, i)), 'color', color);
    end
    hold on;
end

% 设置轴标签
xlabel('$j$', 'interpreter', 'latex');
ylabel('$|\psi_m(j)|^2$', 'interpreter', 'latex');

% 设置图例
lgd = legend([z1,z2,z3],'$\eta=10^{-4}$','$\eta=1$','$\eta=10^{4}$', 'interpreter', 'latex', 'location', 'southeast');
set(lgd, 'Color', [1, 1, 1, 0.5]);

% 设置 x 轴范围
xlim([0, L]);
ylim([10^(-20),1])


end

% 子图5的绘制
function plot_subgraph_5()
    lambda=1.5;
    L=100;
    a=load('psi2_AL_OBC.mat');
    data=a.numericData;
    semilogy(data,'linewidth',3)

    hold on;
    k= 2*log(lambda);

    X=linspace(10,L-10,L-10);
    Y=exp(-k*X);
    Y=Y/sqrt(sum(Y.*Y));

    Y=Y/(10^10);
    semilogy(X,Y,'k--','linewidth',3)

    xlabel('$j$','interpreter','latex')
    ylabel('$|\psi_m(j)|^2$','interpreter','latex')
    lgd=legend('$m=50$','$m=70$','$m=90$','Predicted','interpreter','latex','location','southwest')
    set(lgd, 'Color', [1, 1, 1, 0.5]);
    xlim([0,L])
    ylim([10^(-50),1])
end

% 子图6的绘制
function plot_subgraph_6()
    lambda=1.5;
    L=100;
    a=load('psi2_AL_GBC.mat');
    data=a.numericData;
    semilogy(data,'linewidth',3)
    xlabel('$j$','interpreter','latex')
    ylabel('$|\psi_m(j)|^2$','interpreter','latex')
    lgd=legend('$m=50$','$m=70$','$m=90$','interpreter','latex','location','southeast')
    set(lgd, 'Color', [1, 1, 1, 0.5]);
    xlim([0,L])
    ylim([10^(-50),1])
    
end

% 子图7的绘制
function plot_subgraph_7()
 a=load('GBCEnergyChange1.mat');
 a1=a.Data;
 bc=10.^(-4:0.5:15);
 loglog(bc,a1','linewidth',3)
 
 hold on;
 xlabel('$\eta$','interpreter','latex')
 ylabel('$\overline{\Delta E}$','interpreter','latex')
 xticks([10^-4,10^0,10^5,10^15])
 ylim([10^(-30),10^10])

 % plot the predicted scaling
 
 X=10.^(-3:0.5:5);
 Y=exp(log(X))*10^(-25);
 loglog(X,Y,'k--','linewidth',3)
 
 legend('$L=60$','$L=80$','$L=100$','Predicted','interpreter','latex','location','southeast')
 

end

% 子图8的绘制
function plot_subgraph_8()
 a=load('Eigenstates_AL_bd_gamma1e-04.mat');
Ev=a.Ev;
L=a.L;
color = [0, 0.4470, 0.7410];
for i=5
    semilogy(linspace(1,L,L),Ev(:,i).*conj(Ev(:,i)),'color',color,'linewidth',3)
    hold on;
end

a=load('Eigenstates_AL_bd_gamma1e+00.mat');
Ev=a.Ev;
L=a.L;
color =[0.8500, 0.3250, 0.0980];
for i=5
    semilogy(linspace(1,L,L),Ev(:,i).*conj(Ev(:,i)),'color',color,'linewidth',3)
    hold on;
end

a=load('Eigenstates_AL_bd_gamma1e+04.mat');
Ev=a.Ev;
L=a.L;
color = [0.4660, 0.6740, 0.1880];
for i=5
    semilogy(linspace(1,L,L),Ev(:,i).*conj(Ev(:,i)),'color',color,'linewidth',3)
    hold on;
end

xlabel('$j$','interpreter','latex')
ylabel('$|\psi_m(j)|^2$','interpreter','latex')
lgd=legend('$\eta=10^{-4}$','$\eta=1$','$\eta=10^{4}$','interpreter','latex','location','southeast');
set(lgd, 'Color', [1, 1, 1, 0.5]);
xlim([0,L])


end


function plot_inset1()
%  inset of fig2
    % Adjust the position of the small plot to the center
    a = load('k_with_L.mat');
    Data = a.datak;
    L_all = 50:10:150;
    plot(L_all, Data, '*-', 'linewidth', 2)
    hold on;
    xlabel('$L$', 'interpreter', 'latex')
    ylabel('$k$', 'interpreter', 'latex')
    ylim([9, 9.5])

end


function plot_inset2()
%  inset of fig7
 a=load('GBCEnergyChange1.mat');
 a1=a.Data;
 bc=10.^(-4:0.5:15);
 m=9;
 bc(m)
 y=a1(:,m);
 x=[60,80,100];
 semilogy(x,y,'*-','linewidth',2)
 hold on;

 
 X=60:1:80;
 Y=exp(-log(1.5)*X)*10^(-13);
 semilogy(X,Y,'k--','linewidth',2)
 xticks([60,80,100])
 xlabel('$L$','interpreter','latex')
%  ylabel('$\overline{\Delta E}$','interpreter','latex')
end



