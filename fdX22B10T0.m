% fdX22B10T0 function
% INVERSE HEAT CONDUCTION, Ill-posed problems, 2nd Edition, WILEY, 2023
% Keith A. Woodbury, Hamidreza Najafi, Filippo de Monte, James V. Beck
% INPUTS:
% xd: dimensionless location starting at xd=0 and ending at xd=1
%     xd = x/L
% td: dimensionless time starting at td=0
%     td = alpha*t/L^2
% A:  accuracy desired (1E-A = 10^-A, with A = 2, 3, ..., 15) with respect
%     to the maximum temperature rise that always occurs at xd = 0 and at 
%     the time td assigned for the current X22B10T0 case
% OUPUTS:
% Td: dimensionless temperature calculated at (xd,td) to accuracy A desired
%     Td = (T-Tin)/(q0*L/k)
% Calling Sequence: 
% ierfc(z) = function for computing the complementray error function 
%            integral appearing in the computational analytical solution
function [Td]=fdX22B10T0(xd,td,A)
lengthx=length(xd);
lengtht=length(td);
Td=zeros(lengtht,lengthx); % Preallocating Arrays for speed
for it=1:lengtht  % Begin time loop
    td_it=td(it); % Set current time
    for ix=1:lengthx % Begin space loop 
        xd_ix=xd(ix); % Set current space
        td_dev=(1/(10*A))*(2-xd_ix)^2; % Deviation time
        td_delta=(1/(10*A))*(2+xd_ix)^2; % Second deviation time
        td_qs=(A*log(10)-log(pi^2/2))/pi^2; % Quasi-steady time
        if td_delta > td_qs % This can occurs if A = 2 and xd is close to 1
           td_delta=td_qs;
        end % if td_delta   
        if td_it == 0 % For time t=0 condition
            Td(it,ix)=0; % Set initial temperature
        elseif  td_it <= td_dev 
            % Solution for very small times:
            Td(it,ix)=sqrt(4*td_it)*ierfc(xd_ix/sqrt(4*td_it));
        elseif td_it > td_dev && td_it <= td_delta
            % Solution for small times:
            Td(it,ix)=sqrt(4*td_it)*ierfc(xd_ix/sqrt(4*td_it))+...
                      sqrt(4*td_it)*ierfc((2-xd_ix)/sqrt(4*td_it));
        elseif td_it > td_delta && td_it < td_qs
            % Start solution for large times:
            m_A_max=ceil(sqrt((A*log(10)-log(pi^2/2))/(td_it*pi^2))); % #terms
            Td(it,ix)=td_it+(1/3-xd_ix+xd_ix^2/2); % quasi-steady T solution
            for m=1:m_A_max % Continue solution for large times   
                % Series solutions:
                betam=m*pi; % m-th eigenvalue
                % m-th large-time term of the series solution
                Td(it,ix)=Td(it,ix)-2*exp(-betam^2*td_it)*...
                          cos(betam*xd_ix)/betam^2;
            end % for m
        else
            % Quasi-steady solution for very large times:
            Td(it,ix)=td_it+(1/3-xd_ix+xd_ix^2/2); % quasi-steady temperature solution
        end % if td_it
    end % for ix
end % for it
%
function [ierfc]=ierfc(z)
ierfc =(1/sqrt(pi))*exp(-z^2)-z*erfc(z); % z is any real number