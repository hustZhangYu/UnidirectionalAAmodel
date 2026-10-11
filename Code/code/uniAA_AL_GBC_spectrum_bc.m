clear all;
clc;

% DataE=DataAnalysisA()
% DataAnalysisB()
 DataAnalysisC()

function DataE=DataAnalysisC()

data2=zeros(1,3)

x=10.^(-5:0.5:5);

c=readLegacyCsv('Data1.csv');
semilogx(x,c,'linewidth',4)
hold on;
data2(1,1)=c(11);
c=readLegacyCsv('Data2.csv');
semilogx(x,c,'linewidth',4)
data2(1,2)=c(11);
c=readLegacyCsv('Data3.csv');
semilogx(x,c,'linewidth',4)
data2(1,3)=c(11);

x1=x(1:10);
y1=-60+0.5*linspace(1,10,10)*log(10);
semilogx(x1,y1,'k--','linewidth',2)


xlabel('$\eta$','interpreter','latex','fontsize',24)
ylabel('$\ln(\overline{\Delta E})$','interpreter','latex','fontsize',24)
axis([min(x),max(x),min(c)-25,20])
set(gca,'box','on')
set(gca,'fontsize',18)
    
legend('$L=60$','$L=80$','$L=100$','$k=\ln(10)$','interpreter','latex','location','southeast')


axes('Position',[0.2,0.60,0.2,0.2])
plot([60,80,100],data2,'*-','linewidth',4)
hold on;
set(gca,'fontsize',16)
xlabel('$L$','interpreter','latex','fontsize',16)
ylabel('$\ln(\overline{\Delta E})$','interpreter','latex','fontsize',16)
y1=-35-linspace(1,20,20)*log(1.5);
plot(linspace(61,80,20),y1,'k--','linewidth',2)
end

function DataE=DataAnalysisA()


c=readLegacyCsv('DataERa.csv');
DataE=log(c);
plot(60*ones(1,length(c)),DataE,'.')
hold on;
c=readLegacyCsv('DataERb.csv');
DataE=log(c);
plot(80*ones(1,length(c)),DataE,'.')
c=readLegacyCsv('DataERc.csv');
DataE=log(c);
plot(100*ones(1,length(c)),DataE,'.')
% c=readLegacyCsv('DataERb.csv');
% d=readLegacyCsv('DataEIb.csv');
% DataE(2)=mean(c);
% c=readLegacyCsv('DataERc.csv');
% d=readLegacyCsv('DataEIc.csv');
% DataE(3)=mean(c);


lambda=1.5;
gamma=10^4;
L_all=[60,80,100]
dE=-log(lambda)*L_all+log(gamma)
plot(L_all,dE)

axis([60,100,min(DataE),0])

end



function DataAnalysisB()
a=readLegacyCsv('DataR.csv');
b=readLegacyCsv('DataI.csv');
Data=a+1i*b;

c=readLegacyCsv('DataER.csv');
d=readLegacyCsv('DataEI.csv');
DataE=c+1i*d;

plot(log(c))

L=100;
data=zeros(1,L);
for i=1:L
   psi=Data(:,i); 
   [m,data(i)]=max(log(abs(psi)));
end

m_all=[25,50,75];

for n=1:length(m_all)
    m=m_all(n);
    label=0;
    for k=1:L
       if data(k)==m 
            label=k;
       end   
    end
    psi=Data(:,label);
    plot(log(abs(psi)),'linewidth',4)
    hold on;
end

xlabel('$j$','interpreter','latex','fontsize',24)
ylabel('$\ln(|\psi|)$','interpreter','latex','fontsize',24)
set(gca,'fontsize',18)
legend('$m=25$','$m=50$','$m=75$','interpreter','latex','fontsize',24)
end

function data=readLegacyCsv(filename)
% Recorded CSV inputs live beside the historical C++ prototypes.
dataDir=fullfile(fileparts(mfilename('fullpath')),'..','..','Others','Others');
recordedPath=fullfile(dataDir,filename);
inputPath=fullfile(pwd,filename);
if exist(inputPath,'file')==2
    % Preserve supplied working-directory inputs before recorded examples.
    data=csvread(inputPath,1,0);
else
    data=csvread(recordedPath,1,0);
end
end

