% We use vpa calcualtion to get the eigenstates with different bc
% test

L_all=[60,80,100];
Data=[];
for m=1:length(L_all)
   L=L_all(m);
   data=changeL(L);
   Data=[Data;data]
end

csvwrite( 'GBC.csv','Data');

function Data=changeL(L)
% different boundary condition  
bc=vpa(10.^(-4:0.5:15),50);
                      
Data=zeros(1,length(bc));

for k=1:length(bc)
    omega = vpa((sqrt(5)-1)/2, 50);  % ʹ�� vpa ��߾��ȵ� 50 λ
    lambda = vpa(1.5, 50);  % ʹ�� vpa ��߾��ȵ� 50 λ

    % �߾��Ⱦ�����
    H = diag(ones(1,L-1),-1) + diag(2*lambda*cos(2*pi*omega*linspace(0,L-1,L)));

    % ���� H(1, L) Ԫ��Ϊ 10^(-4)
    H(1,L) = bc(k);  % ʹ�� vpa ȷ������

    % H(1,L) = vpa(0,50); 
    % Ϊ��ȷ������ H �Ǿ�ȷ�ģ�ʹ�� vpa ����ת���ɸ߾�������
    H = vpa(H, 50);  

    % ��������ֵ����������
    [Ev, D] = eig(H);  % ��ʹ�� 'vector'��ֱ�Ӽ�������ֵ����������


    E1=diag(D);

    V=2*lambda*cos(2*pi*omega*linspace(0,L-1,L));


    for m=1:L
        [dE, idx] = min(abs(E1 -V(m)));
        Data(k)=Data(k)+dE;
    end
    
    Data(k)=Data(k)/L;
end




end
