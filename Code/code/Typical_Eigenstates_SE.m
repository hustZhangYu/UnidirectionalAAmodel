% We plot the typical eigenstates in the open boundary conditions.
% The parameter lambda can be changed from 0.5  

L=100;
omega=(sqrt(5)-1)/2;
lambda=0.5;

H=diag(ones(1,L-1),-1)+diag(2*lambda*cos(2*pi*omega*linspace(0,L-1,L)));
H(1,L)=10^(-4);

[Ev,E]=eig(H,'vector');

m_all=[50,70,90];

data=zeros(L,L);
n=1;
for m = 1:L
    semilogy(Ev(:,m).*conj(Ev(:,m)),'linewidth',3)
    hold on;
    data(:,n)=Ev(:,m).*conj(Ev(:,m));
    n=n+1;
end
    
% 将data数组保存为CSV文件
save('psi2_SE_GBC.mat', 'data');

k= 2*log(lambda);

X=linspace(60,L-10,L-10);
Y=exp(-k*X);
Y=Y/sqrt(sum(Y.*Y));

Y=Y/(10^10);
semilogy(X,Y,'k--','linewidth',3)

xlim([0,L])
ylim([10^(-50),1])

xlabel('$j$','interpreter','latex')
ylabel('$|\psi|^2$','interpreter','latex')