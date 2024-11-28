function Qfo=fSFSM(Y, X, r)
%% fSFSM 矩阵最小二乘法
% Ref: Woodbury K A, Najafi H, Monte F D, et al. Inverse Heat Conduction:
% Ill-Posed Problems[M]. 2nd. Hoboken:John Wiley & Sons, 2023.
% Pages 120, Equations 4.66
% Author: Chuansheng LIU
% Data: 20230528
% [Ns, Nt+1] =size(Y)
% [Ns, Nf, Nt]=size(Phi);
% r is is the number of time steps including future data, default is 2.
[Ns, Nf, Nt]=size(X);
Qfo=zeros(Nf, Nt-(r-1));
deltaY=Y(:, 2:end)-Y(:, 1);
dPhi=diff(X, 1, 3);
Xp=[];
for ir=1:r
    Xp=cat(1, Xp, X(:, :, ir));
end
coefficientMatrix=Xp'*Xp;
k=coefficientMatrix^-1*Xp';
for iNt=1:Nt-(r-1)
    That=zeros(Ns, r);
    deltaT=[];
    for ir=1:r
        for j=1:iNt-1
            That(:,ir)=That(:,ir)+dPhi(:, :, iNt-j+ir-1)*Qfo(:, j);
        end
        deltaT=cat(1, deltaT, deltaY(:, iNt+ir-1)-That(:,ir));
    end
    temp=k*deltaT;
    Qfo(:, iNt)=temp;
end
end
