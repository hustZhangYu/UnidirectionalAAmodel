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
    omega = vpa((sqrt(5)-1)/2, 50);  % 使用 vpa 提高精度到 50 位
    lambda = vpa(1.5, 50);  % 使用 vpa 提高精度到 50 位

    % 高精度矩阵构造
    H = diag(ones(1,L-1),-1) + diag(2*lambda*cos(2*pi*omega*linspace(0,L-1,L)));

    % 设置 H(1, L) 元素为 10^(-4)
    H(1,L) = bc(k);  % 使用 vpa 确保精度

    % H(1,L) = vpa(0,50); 
    % 为了确保矩阵 H 是精确的，使用 vpa 函数转换成高精度类型
    H = vpa(H, 50);  

    % 计算特征值和特征向量
    [Ev, D] = eig(H);  % 不使用 'vector'，直接计算特征值和特征向量


    E1=diag(D);

    V=2*lambda*cos(2*pi*omega*linspace(0,L-1,L));


    for m=1:L
        [dE, idx] = min(abs(E1 -V(m)));
        Data(k)=Data(k)+dE;
    end
    
    Data(k)=Data(k)/L;
end




end
