% In this code, we plot all of the eigenstates
%  fixed parameters： disorder strength: 0.5/ boundary condition eta=10^(-4) 
%  differnt system sizes
%  fit ln(|\psi|)=k(j-j0))+c and how k change with L

datak=main();
save('k_with_L.mat','datak')

function datak=main()


L_all=50:10:150;
datak=zeros(1,length(L_all));

for m=1:length(L_all)
   L=L_all(m); 
   omega=(sqrt(5)-1)/2;
   lambda=0.1;

   H=diag(ones(1,L-1),-1)+diag(2*lambda*cos(2*pi*omega*linspace(0,L-1,L)));
   H(1,L)=10^(-4);

   [Ev,E]=eig(H,'vector');
   
   for k=1:L
       psi=Ev(:,k);
       datak(m)=datak(m)+getfit(L,psi); 
   end
   
   datak(m)=datak(m)/L;
end

end


function k0=getfit(L,psi)
% get the wavefunction and return the slope

% 计算波函数的绝对值
psi_abs = abs(psi);

% 选择一个适当的区间进行线性拟合
% 假设我们拟合前半部分数据
x = (1/L:1/L:1)';
y = log(psi_abs(1:L));

% 进行线性拟合，得到斜率
p = polyfit(x, y, 1);
k0 = p(1); % 返回斜率
end
