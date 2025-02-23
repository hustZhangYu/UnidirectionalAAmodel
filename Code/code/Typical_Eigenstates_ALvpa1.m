% AL region 
% generalized boundary condition
% eigenvalue change


L = 100;
omega = vpa((sqrt(5)-1)/2, 50);  % ʹ�� vpa ��߾��ȵ� 50 λ
lambda = vpa(1.5, 50);  % ʹ�� vpa ��߾��ȵ� 50 λ

% �߾��Ⱦ�����
H = diag(ones(1,L-1),-1) + diag(2*lambda*cos(2*pi*omega*linspace(0,L-1,L)));

% ���� H(1, L) Ԫ��Ϊ 10^(-4)
H(1,L) = vpa(0, 50);  % ʹ�� vpa ȷ������

% H(1,L) = vpa(0,50); 
% Ϊ��ȷ������ H �Ǿ�ȷ�ģ�ʹ�� vpa ����ת���ɸ߾�������
H = vpa(H, 50);  

% ��������ֵ����������
[Ev, D] = eig(H);  % ��ʹ�� 'vector'��ֱ�Ӽ�������ֵ����������


E1=diag(double(D));



% ���� m_all ��ֵ
m_all = [50, 70, 90];

% ��ʼ�����ݾ���
data = zeros(L, 3);

Em=[2.4457,-2.9009,1.2715];

n = 1;
for m = m_all
    omgea=(sqrt(5)-1)/2;
    [~, idx] = min(abs(E1 -Em(n) ));
    % �����һ��������������ƽ��ģ
    
    ta = sum(Ev(:,idx).*conj(Ev(:,idx)));
    
    % ʹ�� semilogy ��������������ƽ��ģ
    semilogy(Ev(:,idx).*conj(Ev(:,idx))/ta, 'linewidth', 3)
    hold on;
    
    % �洢����
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

% �� data ת��Ϊ��ֵ����
numericData = double(data);

% ���� data �� CSV �ļ�
save('psi2_AL_OBC.mat', 'numericData');
