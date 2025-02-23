% We plot the typical eigenstates in the open boundary conditions.
% The parameter lambda can be chosen as 1.5

L=100;
omega=(sqrt(5)-1)/2;
lambda=1.5;

H=diag(ones(1,L-1),-1)+diag(2*lambda*cos(2*pi*omega*linspace(0,L-1,L)));
H(1,L)=0;

[Ev,E]=eig(H,'vector');

m_all=[50,70,90];

data=zeros(L,3);
n=1;
for m = m_all
    E(m)
    ta=sum(Ev(:,m).*conj(Ev(:,m)));
    semilogy(Ev(:,m).*conj(Ev(:,m))/ta,'linewidth',3)
    hold on;
    data(:,n)=Ev(:,m).*conj(Ev(:,m));
    n=n+1;
end
    
% 将data数组保存为CSV文件
csvwrite('psi2_AL.csv', 'data');

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