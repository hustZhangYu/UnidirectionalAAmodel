% We plot the typical eigenstates in the open boundary conditions.
% The parameter lambda can be changed from 0.5  

L=100;
omega=(sqrt(5)-1)/2;
lambda=1.5;

H=diag(ones(1,L-1),-1)+diag(2*lambda*cos(2*pi*omega*linspace(0,L-1,L)));
H(1,L)=10^(4);

[Ev,E]=eig(H,'vector');

m_all=[50,70,90];

for m = m_all
    semilogy(Ev(:,m).*conj(Ev(:,m)))
    E(m)
    hold on;
end
    




