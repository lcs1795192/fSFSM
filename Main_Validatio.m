%% This MATLAB script illustrates how to use fSFSM function and proves the same ability as the original file for 1-D IHCP.
% The original MATLAB file is at https://ihc2e.com/wp-content/uploads/codes/chapter4/Example4_2.m
% The detail 1-D IHCP is described in page 123 of Woodbury K A, Najafi H, Monte F D, et al. Inverse Heat Conduction:
% Ill-Posed Problems[M]. 2nd. Hoboken:John Wiley & Sons, 2023.

% Example 4.2 - use finite slab
% 
% Book: INVERSE HEAT CONDUCTION, Ill-posed problems, 2nd Edition, 
% John Wiley & Sons, 2022
% Authors: Keith A. Woodbury, Hamidreza Najafi, 
%          Filippo de Monte, James V. Beck
% INPUTS:
% None.  Various parameters of solution can be changed in this script
% OUPUTS:
% 
% Calling Sequence: 
%  fdX22B20T0
%  fdX22B10T0
clear all
A = 15;         % computational accuracy = 10^(-A)
xs = 0.01;      % meters - depth of the sensor
L = 0.100;       % meters - thickness of the plate
k = 40;         % W/m-K
alpha = 1E-5;   % m^2/s
q_ref = 1;      % W/m^2
qramp = 75000;  % W/m^2/s
Ntimes = 5;    % number of heat flux pulses to estimate
Nround = 3;     % number of decimal places to round Temp 
dt = 5;         % seconds
time = [1:Ntimes]*dt;
% mnon-dimensionalize space and time variables
td = time*alpha/L^2;
xsd = xs/L;
% compute the exact data for this problem
Td = fdX22B20T0( xsd, td, 1*alpha/L^2, A );  % older code
% Td = fdX22B20T0( xsd, td, A )/( 1*alpha/L^2);  % Filippo's code
Temp = [ 30  30+Td'*qramp*L/k ];
Temp = round(Temp,Nround);
Y = Temp(2:length(Temp));  % only the last are needed
Y0 = Temp(1);   % this is the initial condition temperature
%
% now non-dimensionalize temperature data 
Yd = (Y - Y(1))/(q_ref*L/k);
% compute the responses for the X22B10T0 problem
Tresp = fdX22B10T0( xsd, td, A);
% these are the dimensionsless phis
phid = Tresp;
% compute the dphis
dphid =  [ Tresp(2:length(time))- Tresp(1:length(time)-1) ];
% add dimensions
phi = phid * L/k;
dphi = dphid * L/k;
%
%  now use FS to solve
%
rvec=[2:2];
Ydata=Temp(2:Ntimes+1)-Temp(1);

for r=rvec
    clear K0
    K0 = zeros( r, 1);
    Ksum = sum( phi(1:r).^2 );
    K0 = phi(1:r)/Ksum;
    qhat_FS = zeros(1,Ntimes-r+1);
    for M = 1:Ntimes-r+1
        That_0 = zeros(r,1);    % temporary variable
        for ir = 1:r
            for k = 1:M-1
                That_0(ir) = That_0(ir) + qhat_FS(k)*dphi(M+ir-k-1);
%                fprintf(' M=%d  ir=%d  M+ir-k-1=%d \n',M,ir,M+ir-k-1);
            end
            qhat_FS(M) = qhat_FS(M) + ...
                               (Ydata(M+ir-1) - That_0(ir))*K0(ir);
        end
    end
end
% compute exact solution at time -dt/2 for best comparison
time-time(1)/2;
qexact = qramp*(time-time(1)/2);

%% This line is the only line different form the original file, which shows the usage of fSFSM.
% Qfo shoule be the same as qhat_FS
Qfo=fSFSM(Temp,reshape(phi,1,1,5),2);
%% End of modifications.

clear err_FS;
err_FS = (1 - qhat_FS./qexact(1:Ntimes-r+1))*100;
fprintf(' Temps (C):\n');
fprintf(' %7.3f ',Temp);
fprintf('\n')
fprintf(' Solution by Function Specification with piecewise constant q\n')
fprintf(' K (W/m^2/C):\n');
fprintf(' %6.4e ',K0);
fprintf('\n')
fprintf(' qhat_FS (W/m^2):\n');
fprintf(' %6.1f ',qhat_FS);
fprintf('\n')
fprintf(' qexact (W/m^2):\n');
fprintf(' %6.1f ',qexact);
fprintf('\n')
fprintf(' relative errors:\n');
fprintf(' %7.1f%% ',err_FS);
fprintf('\n')
