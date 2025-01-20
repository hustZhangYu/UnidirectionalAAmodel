
Main_Plot()

function H=Main_Plot()
gamma_all=[0.0001,1,10000];
% 定义不同的标记符号
markers = ['o', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '*']; 
Lines = {'-', '--', '-.'};
n=1;
for gamma= gamma_all
    marker=markers(n);
    line=Lines{n};


    n=n+1;
    Plot_SFL_Spectrum(gamma,marker,line)
    H=Plot_AL_Spectrum(gamma,marker);
end

legend('Analytical Results $\lambda=0.5,\eta=10^{-4}$','$\lambda=0.5,\eta=10^{-4}$','$\lambda=1.5,\eta=10^{-4}$','Analytical Results $\lambda=0.5,\eta=10^{0}$','$\lambda=0.5,\eta=10^{0}$','$\lambda=1.5,\eta=10^{0}$','Analytical Results $\lambda=0.5,\eta=10^{4}$','$\lambda=0.5,\eta=10^{4}$','$\lambda=1.5,\eta=10^{4}$','interpreter','latex')
end


function H=Plot_AL_Spectrum(gamma1,marker1)

lambda1=1.5;

L=100;
lambda=lambda1;

H=Ham(L,lambda,gamma1);
[Ev,E]=eig(H,'vector');
hold on;
a1=plot(imag(E),real(E),marker1,'markersize',5,'linewidth',2);

filename = sprintf('AL_Spectrum_gamma_%.4f.mat', gamma1);
save(filename, 'E');

xlabel('$\Im(E)$','interpreter','latex','Fontsize',24)
ylabel('$\Re(E)$','interpreter','latex','Fontsize',24)
view(90,90)

hold on

end


function Plot_SFL_Spectrum(gamma1,marker1,line)

lambda1=0.5;

L=100;


 lambda=lambda1;
% 先给出解析求解李雅普诺夫指数的结果图
a=-3.5:0.1:3.5;
b=-3:0.1:3;
data=zeros(length(a),length(b));
for m=1:length(a)
    for n=1:length(b)
          a1=a(m);b1=b(n);
          c=a1+1i*b1;
          r1= a1^2/(1+(2*lambda)^2/(4*exp(2*log(gamma1)/L)))^2+b1^2/(1-(2*lambda)^2/(4*exp(2*log(gamma1)/L)))^2;
          data(m,n)=r1;
    end
end

contour(b,a,real(data),[exp(2*log(gamma1)/L),exp(2*log(gamma1)/L)],line,'linewidth',2,'color','k')

lambda=lambda1;

H=Ham(L,lambda,gamma1);
[Ev,E]=eig(H,'vector');
hold on;
a1=plot(imag(E),real(E),marker1,'markersize',5,'linewidth',2);

% 保存 data 数据为 CSV 文件
filename = sprintf('SFL_Spectrum_gamma_%.4f.mat', gamma1);
filename1= 'SFL_Spectrum_Analytic.mat';
save(filename, 'E');
save(filename1, 'data');

xlabel('$\Im(E)$','interpreter','latex','Fontsize',24)
ylabel('$\Re(E)$','interpreter','latex','Fontsize',24)
view(90,90)

hold on

end


function H=Ham(L,lambda,gamma)
% 我们给出相应的哈密顿量
    omega=(sqrt(5)-1)/2;
    V=2*lambda*diag(cos(2*pi*omega*linspace(1,L,L)));
    H=diag(ones(1,L-1),-1)+V;
    H(1,L)=gamma;
end