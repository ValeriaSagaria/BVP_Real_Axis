fprintf('script vecchio metodo nuovo\n');
h=@(x) (abs(x).^(13/2).*exp(-3*x.^2/2));
kind=2;
tetaOld=0.5;
tetaNy=tetaOld;
ANew=@(x) (1./(x.^8+5).*exp(x.^2/2));  
t=linspace(-6,6,30000);
t=[-6 2 6];
t=6;
de=1/2+0.2;
%  t=2;
%% Metodi per BVP
l = 9;
MNew=zeros(l,length(t)+3);
for i=2:l
    mm=2^i;
%     [m,jNew,CondNew,FmNew] = BVPNew(mm,h,ANew,t,kind);
    [FmNew,CondNew,jNew] = NystromMethodFIE(mm,kind,h,ANew,t);
    MNew(i-1,:)=[m jNew CondNew FmNew];
end

%% Soluzione esatta
mm=600;
% [m,jNew,CondNew,FmNew] = BVPNew(mm,h,ANew,t,kind);
 [FmNew,CondNew,jNew] = NystromMethodFIE(mm,kind,h,ANew,t);
MNew(l,:)=[m jNew CondNew FmNew]; % Esatta
% fprintf('Condizionamento e troncamento soluzione esatta nuovo metodo \n');
% CondNew
% jNew

MNew(:,2:end)

%% Errore assoluto
merrNew=zeros(l-1,2);
errNew=zeros(length(t),l-1);

for k=2:l
    for i=1:length(t)
        errNew(i,k-1)=abs(MNew(l,i+3)-MNew(k-1,i+3));
    end
    merrNew(k-1,:)=[2^k max(errNew(:,k-1))];
end

%% Tabelle m & j & cond(A_m) & err_m
% fprintf('Errori nuovo metodo \n');
% [MNew(1:l-1,1:3) merrNew(:,2)]