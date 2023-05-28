function Qfo=fSFSM(Y, Phi, r)
%% fSFSM 矩阵最小二乘法
% Ref: Woodbury K A, Najafi H, Monte F D, et al. Inverse Heat Conduction:
% Ill-Posed Problems[M]. 2nd. Hoboken:John Wiley & Sons, 2023.
% Pages 120, Equations 4.66
% Author: Chuansheng LIU
% Data: 20230528
% [Ns, Nt+1] =size(Y)
% [Ns, Nf, Nt]=size(Phi);
% r is a scalar.
[Ns, Nf, Nt]=size(Phi);
Qfo=zeros(Nf, Nt-(r-1));
deltaY=Y(:, 2:end)-Y(:, 1);
dPhi=diff(Phi, 1, 3);
X=[];
for ir=1:r
    X=cat(1, X, Phi(:, :, ir));
end
coefficientMatrix=X'*X;
k=coefficientMatrix^-1*X';
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
