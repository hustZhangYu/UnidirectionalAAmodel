EigenstatesPlot1()



function EigenstatesPlot1()

digits(50); % 设置变量精度为 50 位

L = 100;

E0=-2.9009;

E0=2.4457;

% E0=1.2715;
% 使用 vpa 定义 gamma_all
gamma_all = vpa([10^(-4), 10^(0), 10^(4)]);

for m = 1:length(gamma_all)
    gamma1 = gamma_all(m); % 使用 vpa 的 gamma1
    lambda = vpa(1.5);     % 将 lambda 也设置为高精度
    
    % 假设 Ham 是用户定义的函数，返回高精度的结果
    H = vpa(Ham(L, lambda, gamma1));
    [Ev, E] = eig(H); % Ev 和 E 都是高精度计算结果
    
    [r, s] = GetR(L, Ev); % 假设 GetR 函数返回所需的高精度结果
    
    log_gamma1 = log(gamma1) / L; % 高精度计算
    
    num_curves = L;
    if m == 1
        color = [0, 0.4470, 0.7410]; % 蓝色
    elseif m == 2
        color = [0.8500, 0.3250, 0.0980]; % 红色
    elseif m == 3
        color = [0.4660, 0.6740, 0.1880]; % 绿色
    end
    
    alphas = linspace(0.1, 1, num_curves); % 透明度从 0.1 到 1
    
    Ed=double(diag(E));
    [~,idx]=min(abs(Ed-E0));
    for i = idx
        psi = Ev(:, i).*conj(Ev(:, i)); % 高精度计算
        if m == 1
            a = semilogy(linspace(1, L, L) / L, psi, 'Color', color, 'marker', 'o', 'linewidth', 2);
        elseif m == 2
            b = semilogy(linspace(1, L, L) / L, psi, 'Color', color, 'marker', 's', 'linewidth', 2);
        elseif m == 3
            c = semilogy(linspace(1, L, L) / L, psi, 'Color', color, 'marker', 'p', 'linewidth', 2);
        end
        hold on;
    end
    
    % 保存结果，文件名保持一致
    filename = sprintf('Eigenstates_AL_bd_gamma%.0e.mat', double(gamma1)); % 转换为双精度以保存
    save(filename, 'Ev', 'gamma1', 'lambda', 'L');
    
    filename = sprintf('Eigenvalues_AL_bd_gamma%.0e.mat', double(gamma1)); % 转换为双精度以保存
    save(filename, 'Ed', 'gamma1', 'lambda', 'L');
    
    
    xlabel('$j/L$', 'interpreter', 'latex');
    ylabel('$|\psi|^2$', 'interpreter', 'latex');
    
    set(gca, 'fontsize', 18);
end

% 图例
legend([a, b, c], '$\eta=10^{-4}$', '$\eta=1$', '$\eta=10^{4}$', 'interpreter', 'latex');



end



function EigenstatesPlot()

L=100;

gamma_all=[10^(-4),10^(0),10^(4)];

for m=1:length(gamma_all)
    gamma1=gamma_all(m);
    %  gamma1=1;
    %  gamma1=10^4;
    lambda=1.5;    
    
    H=Ham(L,lambda,gamma1);
    [Ev,E]=eig(H,'vector');
    
    [r,s]=GetR(L,Ev);
    
    log(gamma1)/L;
    
    num_curves = L;
    if m==1
    color = [0, 0.4470, 0.7410]; % 蓝色
    end
    if m==2
    color =[0.8500, 0.3250, 0.0980]; % 红色
    end
    if m==3
    color = [0.4660, 0.6740, 0.1880]; %  绿色
    end
    alphas = linspace(0.1, 1, num_curves); % 透明度从0.1到1
    
    for i=5
        psi=log(abs(Ev(:,i)));
        if m==1
        a=plot(linspace(1,L,L)/L,psi,'Color', color,'marker','o','linewidth',2);
        end
        if m==2
        b=plot(linspace(1,L,L)/L,psi,'Color', color,'marker','s','linewidth',2);
        end
        if m==3
        c=plot(linspace(1,L,L)/L,psi,'Color', color,'marker','p','linewidth',2);
        end
        hold on;
    end
    
    filename = sprintf('Eigenstates_AL_bd_gamma%.0e.mat', gamma1);
    save(filename, 'Ev', 'gamma1', 'lambda', 'L');
    
    xlabel('$j/L$','interpreter','latex')
    ylabel('$\ln(|\psi|)$','interpreter','latex')
    
    
    set(gca,'fontsize',18)
end


legend([a,b,c],'$\eta=10^{-4}$','$\eta=1$','$\eta=10^{4}$','interpreter','latex')

end

function [r,s]=GetR(L,Ev)

SlopeData=zeros(L,2);

Xdata=[];
Ydata=[];

for m=1:L
    psi=Ev(:,m);
    y1=log(abs(psi));
    y=y1';
    a1=linspace(1,L,L)/L;

    data11=isinf(y);
    [a ]=find(data11==1);
    y(a)=[];
    a1(a)=[];
   
    data11=isnan(y);
    [a ]=find(data11==1);
    y(a)=[];
    a1(a)=[];
    
    Xdata=[Xdata,a1];
    Ydata=[Ydata,y];

%     SlopeData(m,1)=r;
%     SlopeData(m,2)=s;

end

[r,s]=Slope(Xdata,Ydata);
% r=mean(SlopeData(:,1));
% s=mean(SlopeData(:,2));
 
end

function [a,b]=Slope(x,y)
% This code is used to take linerfitting 
% input the x, y data and return the slope


r=polyfit(x,y,1);
a=r(1);
b=r(2);

end

function testEnergyPlot()


% 先给出解析求解李雅普诺夫指数的结果图
a=-3.5:0.1:3.5;
b=-3:0.1:3;
data=zeros(length(a),length(b));
for m=1:length(a)
    for n=1:length(b)
          a1=a(m);b1=b(n);
          c=a1+1i*b1;
          if a1>0
              data(m,n)=log((c+sqrt(c*c-1))/2);
          end
          if a1<0
              data(m,n)=log((c-sqrt(c*c-1))/2);
          end
          
          if a1==0
              if b1>0
                  data(m,n)=log((c+sqrt(c*c-1))/2);
              end
              if b1<=0
                   data(m,n)=log((c-sqrt(c*c-1))/2);
              end
          end
    end
end

gamma1=0.001;
gamma2=10;
gamma3=1000;
L=100;
lambda=0.5;

 figure()
contour(b,a,real(data),[log(gamma1)/L-log(2*lambda),log(gamma1)/L-log(2*0.2)],'linewidth',2,'color','k')
hold on;
contour(b,a,real(data),[log(gamma2)/L-log(2*lambda),log(gamma2)/L-log(2*0.2)],'linewidth',2,'color','k')
contour(b,a,real(data),[log(gamma3)/L-log(2*lambda),log(gamma3)/L-log(2*0.2)],'linewidth',2,'color','k')



H=Ham(L,lambda,gamma1);
[Ev,E]=eig(H,'vector')
hold on;
a1=plot(imag(E)/(2*lambda),real(E)/(2*lambda),'.','markersize',15)

H=Ham(L,lambda,gamma2);
[Ev,E]=eig(H,'vector')
hold on;
a2=plot(imag(E)/(2*lambda),real(E)/(2*lambda),'.','markersize',15)

H=Ham(L,lambda,gamma3);
[Ev,E]=eig(H,'vector')
hold on;
a3=plot(imag(E)/(2*lambda),real(E)/(2*lambda),'.','markersize',15)

lambda=0.2;
H=Ham(L,lambda,gamma1);
[Ev,E]=eig(H,'vector')
hold on;
a4=plot(imag(E)/(2*lambda),real(E)/(2*lambda),'.','markersize',15)

H=Ham(L,lambda,gamma2);
[Ev,E]=eig(H,'vector')
hold on;
a5=plot(imag(E)/(2*lambda),real(E)/(2*lambda),'.','markersize',15)

H=Ham(L,lambda,gamma3);
[Ev,E]=eig(H,'vector')
hold on;
a6=plot(imag(E)/(2*lambda),real(E)/(2*lambda),'.','markersize',15)

xlabel('$\Im(E)$','interpreter','latex','Fontsize',24)
ylabel('$\Re(E)$','interpreter','latex','Fontsize',24)
view(90,90)
legend([a1,a2,a3,a4,a5,a6],'$\lambda=0.5,\gamma=10^{-4}$','$\lambda=0.5,\gamma=1$','$\lambda=0.5,\gamma=10^{4}$','$\lambda=0.2,\gamma=10^{-4}$','$\lambda=0.2,\gamma=1$','$\lambda=0.2,\gamma=10^{4}$','interpreter','latex','Fontsize',24,'location','northeast')

end

function testEnergy()


% 先给出解析求解李雅普诺夫指数的结果图
a=-5.5:0.1:5.5;
b=-5:0.1:5;
data=zeros(length(a),length(b));
for m=1:length(a)
    for n=1:length(b)
          a1=a(m);b1=b(n);
          c=a1+1i*b1;
          if a1>0
              data(m,n)=log((c+sqrt(c*c-1))/2);
          end
          if a1<0
              data(m,n)=log((c-sqrt(c*c-1))/2);
          end
          
          if a1==0
              if b1>0
                  data(m,n)=log((c+sqrt(c*c-1))/2);
              end
              if b1<=0
                   data(m,n)=log((c-sqrt(c*c-1))/2);
              end
          end
    end
end

gamma1=0.1;
gamma2=10;
gamma3=100;
L=100;
lambda=0.8;

 figure()
contour(b,a,real(data),[log(gamma1)/L-log(2*lambda),log(gamma1)/L-log(2*0.2),log(gamma1)/L-log(2*0.1)])
hold on;
contour(b,a,real(data),[log(gamma2)/L-log(2*lambda),log(gamma2)/L-log(2*0.2),log(gamma2)/L-log(2*0.1)])
contour(b,a,real(data),[log(gamma3)/L-log(2*lambda),log(gamma3)/L-log(2*0.2),log(gamma3)/L-log(2*0.1)])
xlabel('imag')
ylabel('real')
view(90,90)


H=Ham(L,lambda,gamma1);
[Ev,E]=eig(H,'vector')
hold on;
plot(imag(E)/(2*lambda),real(E)/(2*lambda),'.')

H=Ham(L,lambda,gamma2);
[Ev,E]=eig(H,'vector')
hold on;
plot(imag(E)/(2*lambda),real(E)/(2*lambda),'.')

H=Ham(L,lambda,gamma3);
[Ev,E]=eig(H,'vector')
hold on;
plot(imag(E)/(2*lambda),real(E)/(2*lambda),'.')

lambda=0.2;
H=Ham(L,lambda,gamma1);
[Ev,E]=eig(H,'vector')
hold on;
plot(imag(E)/(2*lambda),real(E)/(2*lambda),'.')

H=Ham(L,lambda,gamma2);
[Ev,E]=eig(H,'vector')
hold on;
plot(imag(E)/(2*lambda),real(E)/(2*lambda),'.')

H=Ham(L,lambda,gamma3);
[Ev,E]=eig(H,'vector')
hold on;
plot(imag(E)/(2*lambda),real(E)/(2*lambda),'.')


lambda=0.1;
H=Ham(L,lambda,gamma1);
[Ev,E]=eig(H,'vector')
hold on;
plot(imag(E)/(2*lambda),real(E)/(2*lambda),'.')

H=Ham(L,lambda,gamma2);
[Ev,E]=eig(H,'vector')
hold on;
plot(imag(E)/(2*lambda),real(E)/(2*lambda),'.')

H=Ham(L,lambda,gamma3);
[Ev,E]=eig(H,'vector')
hold on;
plot(imag(E)/(2*lambda),real(E)/(2*lambda),'.')

end

function Eigenvalues()

L=100;
lambda=0.9;
gamma=0.5;
H=Ham(L,lambda,gamma);
[Ev,E]=eig(H,'vector')

plot(real(E)/(2*lambda),imag(E)/(2*lambda),'.')
axis([-3.5,3.5,-3,3])
hold on;

end

function Eigenstates()

L=100;
lambda=0.9;
gamma=0.1;
H=Ham(L,lambda,gamma);
[Ev,E]=eig(H,'vector')

m=10;

psi=Ev(:,m);

plot(log10(psi.*conj(psi)))

end



function H=Ham(L,lambda,gamma)
% 我们给出相应的哈密顿量
    omega=(sqrt(5)-1)/2;
    V=2*lambda*diag(cos(2*pi*omega*linspace(0,L-1,L)));
    H=diag(ones(1,L-1),-1)+V;
    H(1,L)=gamma;
end