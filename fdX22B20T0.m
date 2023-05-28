% fdX22B20T0 function
% INVERSE HEAT CONDUCTION, Ill-posed problems, 2nd Edition, WILEY, 2023
% Keith A. Woodbury, Hamidreza Najafi, Filippo de Monte, James V. Beck
% INPUTS:
% xd: dimensionless location starting at xd=0 and ending at xd=1
%     xd = x/L
% td: dimensionless time starting at td=0
%     td = alpha*t/L^2
% trefd: dimensionless reference time
%     trefd = alpha*tref/L^2
% A:  accuracy desired (1E-A = 10^-A, with A = 2, 3, ..., 15) with respect
%     to the maximum temperature rise that always occurs at xd = 0 and at 
%     the time td assigned for the current X22B20T0 case
% OUPUTS:
% Td: dimensionless temperature calculated at (xd,td) to accuracy A desired
%     Td = (T-Tin)/(q0*L/k)
% Calling Sequence: 
% i3erfc(z) = function for computing the complementray error function 
%             integral appearing in the computational analytical solution
function [Td]=fdX22B20T0(xd,td,trefd,A)
lengthx=length(xd);
lengtht=length(td);
Td=zeros(lengtht,lengthx); % Preallocating Arrays for speed
for it=1:lengtht  % Begin time loop
    td_it=td(it); % Set current time
    for ix=1:lengthx % Begin space loop 
        xd_ix=xd(ix); % Set current space
            td_d=(1/(10*A))*(2-xd_ix)^2; % Deviation time
        td_delta=(1/(10*A))*(2+xd_ix)^2; % Second deviation time
        td_qs=(A*log(10)-log(pi^4/2))/pi^2; % Quasi-steady time
        if td_delta > td_qs % This occurs only if A = 2 and concerns any xd 
           td_delta=td_qs;
        end % if td_delta   
        if td_it == 0 % For time t=0 condition
            Td(it,ix)=0; % Set initial temperature
        elseif  td_it <= td_d 
            % Solution for very small times:
            Td(it,ix)=(4*td_it)^(3/2)*i3erfc(xd_ix/sqrt(4*td_it))/trefd;
        elseif td_it > td_d && td_it <= td_delta
            % Solution for small times:
            Td(it,ix)=((4*td_it)^(3/2)*i3erfc(xd_ix/sqrt(4*td_it))+...
                      (4*td_it)^(3/2)*i3erfc((2-xd_ix)/sqrt(4*td_it)))/...
                      trefd;
        elseif td_it > td_delta && td_it < td_qs
            % Start solution for large times
            m_A_max=ceil(sqrt((A*log(10)-log(pi^4/2))/(td_it*pi^2))); % #terms
            Td(it,ix)=((td_it^2/2)+(1/3-xd_ix+xd_ix^2/2)*td_it+...
                      (-1/45+xd_ix^2/6-xd_ix^3/6+xd_ix^4/24))/...
                      trefd; % quasi-steady temperature solution
            for m=1:m_A_max % Continue solution for large times   
                betam=m*pi; % m-th eigenvalue
                % m-th large-time term of the series solution
                Td(it,ix)=Td(it,ix)+(2/trefd)*exp(-betam^2*td_it)*...
                          cos(betam*xd_ix)/betam^4;
            end % for m
        else
            % Quasi-steady solution for very large times:
            Td(it,ix)=((td_it^2/2)+(1/3-xd_ix+xd_ix^2/2)*td_it+...
                      (-1/45+xd_ix^2/6-xd_ix^3/6+xd_ix^4/24))/...
                      trefd; % quasi-steady temperature solution
        end % if td_it
    end % for ix
end % for it
%
function [i3erfc]=i3erfc(z)
i3erfc =1/12*(2*exp(-z^2)*(1+z^2)/sqrt(pi)-erfc(z)*z*(3+2*z^2)); % z is any real number