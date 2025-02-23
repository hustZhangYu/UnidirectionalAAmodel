figure()
colors = {[1 0 0], [0 1 0], [1 0.5 0], [1 0 1]};

b=load('eigenvalues.mat');
E=b.E;
subplot(1,2,1)
plot(real(E),imag(E),'.','markersize',20)
hold on;
xlim([-4,4])
ylim([-2,2])


a=load('eigenvectors.mat');
Data=a.Ev;
Data=Data';

z=[4,7,15,20];
[E1,index]=sort(real(E));

marker_all=['o','s','d','p'];

for i=1:4
    m=z(i);
   psi=Data(:,index(m));
   subplot(1,2,2)
   semilogy(abs(psi),'linewidth',2,'color',colors{i},'marker',marker_all(i));
    hold on;
    subplot(1,2,1)
    plot(real(E(index(m))),imag(E(index(m))),'o','markersize',8,'color',colors{i},'marker',marker_all(i),'linewidth',2)
    hold on;
end

subplot(1,2,1)
xlabel('$\Re(E)$','interpreter','latex','fontsize',18)
ylabel('$\Im(E)$','interpreter','latex','fontsize',18)
set(gca,'fontsize',18)

lambda=2;
eta=10^(-4);
% y=-2:0.1:2;
kappa=2;

L=40;

y=-2:0.1:2;
x=((lambda/eta^(2/L))^(1/(kappa-1)))*ones(1,length(y));
x=ones(1,length(y))./x;


eta^(kappa/(L*(kappa-1)))/lambda^(kappa-1)

a=plot(x,y,'k--','linewidth',3)
plot(-x,-y,'k--','linewidth',3)
legend(a,'ME')
%     plot(-x,y,'b--','linewidth',2)




subplot(1,2,2)
xlabel('$j$','interpreter','latex','fontsize',18)
ylabel('$\ln(|\psi_j|)$','interpreter','latex','fontsize',18)
set(gca,'fontsize',18)