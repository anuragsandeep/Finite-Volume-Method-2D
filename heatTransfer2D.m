%% Heat Transfer in 2-D using line by line visitation
%% with Gauss-Seidel and TDMA 
%% Anurag Sandeep K. 

%clear all 
%clc
function [T,count] = heatTransfer2D(NX,NY,alpha)

% INPUT PARAMETERS
xlen=1; % x-dimension
ylen=1; % y-dimension
k=386;  % heat transfer co-efficient [W/m/K]
J=10000; % heat flux
ITMAX=NX;
JTMAX=NY;
dx=xlen/(ITMAX-2); % grid-size in x
dy=ylen/(JTMAX-2); % grid-size in y
T=zeros(ITMAX,JTMAX);

% BOUNDARY CONDITIONS
T(:,1)=50;
T(1,:)=100;
T(end,:)=50;
T(:,end)=-J*(dx/(2*k)) + T(:,end-1);
A=[];B=[];
eps=1e-5;
error=eps+1;
count=0;

while (error > eps)
    count=count+1;

%--------------------------------------------------------------------------
% L E F T - R I G H T (LEFT COLUMN)  ------->
%--------------------------------------------------------------------------
% 1. Equations for bottom-left CV
A(1,2)=(k*(3*(dy/dx)+3*(dx/dy)))/alpha;
A(1,3)=-k*dx/dy;
B(1)=((1-alpha)/alpha)*k*(3*(dy/dx)+3*(dx/dy))*T(end-1,2)+...
    2*k*(dy/dx)*T(end-1,1)+2*k*(dy/dx)*T(end,2)+k*(dy/dx)*T(end-1,3);

% 2. Equations for left-most internal CV's
for jj=3:JTMAX-2
    A(jj-1,1)=-k*(dx/dy);
    A(jj-1,2)=(k*(3*(dy/dx)+2*(dx/dy)))/alpha;
    A(jj-1,3)=-k*(dx/dy);
    B(jj-1)=((1-alpha)/alpha)*k*((3*(dy/dx)+2*(dx/dy))*T(end-jj+1,2))+...
    2*k*(dy/dx)*T(end-jj+1,1)+k*(dy/dx)*T(end-jj+1,3);
end

% 3. Equations for top-left CV
A(JTMAX-2,1)=-k*dx/dy;
A(JTMAX-2,2)=(k*(3*(dy/dx)+3*(dx/dy)))/alpha;
B(JTMAX-2)=((1-alpha)/alpha)*k*(3*(dy/dx)+3*(dx/dy))*T(2,2)+...
    2*k*(dy/dx)*T(2,1)+2*k*(dy/dx)*T(1,2)+k*(dy/dx)*T(2,3);

T(end-1:-1:2,2)=Tridiagonal(JTMAX-2,A,B);

% INTERNAL LEFT-RIGHT SWEEP
for ix=3:ITMAX-2
    % 1. Equations for bottom CV
    A(1,2)=(k*(2*(dy/dx)+3*(dx/dy)))/alpha;
    A(1,3)=-k*dx/dy;
    B(1)=((1-alpha)/alpha)*k*(2*(dy/dx)+3*(dx/dy))*T(end-1,ix)+...
    k*(dy/dx)*T(end-1,ix-1)+2*k*(dy/dx)*T(end,ix)+k*(dy/dx)*T(end-1,ix+1);

    % 2. Equations for left-most internal CV's
    for jj=3:JTMAX-2
        A(jj-1,1)=-k*(dx/dy);
        A(jj-1,2)=(k*(2*(dy/dx)+2*(dx/dy)))/alpha;
        A(jj-1,3)=-k*(dx/dy);
        B(jj-1)=((1-alpha)/alpha)*k*((2*(dy/dx)+2*(dx/dy))*T(end-jj+1,ix))+...
        k*(dy/dx)*T(end-jj+1,ix-1)+k*(dy/dx)*T(end-jj+1,ix+1);
    end

    % 3. Equations for top-left CV
    A(JTMAX-2,1)=-k*dx/dy;
    A(JTMAX-2,2)=(k*(2*(dy/dx)+3*(dx/dy)))/alpha;
    B(JTMAX-2)=((1-alpha)/alpha)*k*(2*(dy/dx)+3*(dx/dy))*T(2,ix)+...
    k*(dy/dx)*T(2,ix-1)+2*k*(dy/dx)*T(1,ix)+k*(dy/dx)*T(2,ix+1);

    T(end-1:-1:2,ix)=Tridiagonal(JTMAX-2,A,B);    
end

% L E F T - R I G H T (RIGHT COLUMN)
% 1. Equations for bottom-left CV
A(1,2)=(k*(3*(dy/dx)+3*(dx/dy)))/alpha;
A(1,3)=-k*dx/dy;
B(1)=((1-alpha)/alpha)*k*(3*(dy/dx)+3*(dx/dy))*T(end-1,JTMAX-1)+...
    k*(dy/dx)*T(end-1,JTMAX-2)+2*k*(dy/dx)*T(end,JTMAX-1)+2*k*(dy/dx)*T(end-1,JTMAX);

% 2. Equations for right-most internal CV's
for jj=3:JTMAX-2
    A(jj-1,1)=-k*(dx/dy);
    A(jj-1,2)=(k*(2*(dy/dx)+3*(dx/dy)))/alpha;
    A(jj-1,3)=-k*(dx/dy);
    B(jj-1)=((1-alpha)/alpha)*k*((3*(dy/dx)+2*(dx/dy))*T(end-jj+1,JTMAX-1))+...
    k*(dy/dx)*T(end-jj+1,JTMAX-2)+2*k*(dy/dx)*T(end-jj+1,JTMAX);
end

% 3. Equations for top-left CV
A(JTMAX-2,1)=-k*dx/dy;
A(JTMAX-2,2)=(k*(3*(dy/dx)+3*(dx/dy)))/alpha;
B(JTMAX-2)=((1-alpha)/alpha)*k*(3*(dy/dx)+3*(dx/dy))*T(2,JTMAX-1)+...
    k*(dy/dx)*T(2,JTMAX-2)+2*k*(dy/dx)*T(1,JTMAX-1)+2*k*(dy/dx)*T(2,JTMAX);

T(end-1:-1:2,JTMAX-1)=Tridiagonal(JTMAX-2,A,B);

% Update the right boundary
T(:,end)=-J*(dx/(2*k)) + T(:,end-1);


%--------------------------------------------------------------------------
% S O U T H - N O R T H (BOTTOM ROW)    ^  
%                                       |
%                                       |
%--------------------------------------------------------------------------
% BOTTOM-MOST 
% 1. Equation for the bottom left-most CV
A(1,2)=(k*(3*(dy/dx)+3*(dx/dy)))/alpha;
A(1,3)=-k*dy/dx;
B(1)=((1-alpha)/alpha)*k*(3*(dy/dx)+3*(dx/dy))*T(end-1,2)+...
    2*k*(dy/dx)*T(end-1,1)+2*k*(dy/dx)*T(end,2)+k*(dy/dx)*T(end-2,2);

% 2. Equations for bottom-most internal CV's
for ii=3:ITMAX-2
    A(ii-1,1)=-k*(dx/dy);
    A(ii-1,2)=(k*(2*(dy/dx)+3*(dx/dy)))/alpha;
    A(ii-1,3)=-k*(dx/dy);
    B(ii-1)=((1-alpha)/alpha)*k*((3*(dy/dx)+2*(dx/dy))*T(end-1,ii))+...
    2*k*(dy/dx)*T(end,ii)+k*(dy/dx)*T(end-2,ii);
end

% 3. Equations for bottom right-most CV
A(ITMAX-2,1)=-k*dx/dy;
A(ITMAX-2,2)=(k*(3*(dy/dx)+3*(dx/dy)))/alpha;
B(ITMAX-2)=((1-alpha)/alpha)*k*(3*(dy/dx)+3*(dx/dy))*T(end-1,end-1)+...
    k*(dy/dx)*T(end-2,end-1)+2*k*(dy/dx)*T(end,end-1)+2*k*(dy/dx)*T(end-1,end);

T(end-1,2:end-1)=Tridiagonal(ITMAX-2,A,B);

% INTERNAL SOUTH-NORTH SWEEP
for jy=3:JTMAX-2
    % 1. Equation for the bottom left-most CV
    A(1,2)=(k*(3*(dy/dx)+2*(dx/dy)))/alpha;
    A(1,3)=-k*dy/dx;
    B(1)=((1-alpha)/alpha)*k*(3*(dy/dx)+2*(dx/dy))*T(end-jy+1,2)+...
    2*k*(dy/dx)*T(end-jy+1,1)+k*(dy/dx)*T(end-jy+2,2)+k*(dy/dx)*T(end-jy,2);

    % 2. Equations for bottom-most internal CV's
    for ii=3:JTMAX-2
    A(ii-1,1)=-k*(dx/dy);
    A(ii-1,2)=(k*(2*(dy/dx)+2*(dx/dy)))/alpha;
    A(ii-1,3)=-k*(dx/dy);
    B(ii-1)=((1-alpha)/alpha)*k*((2*(dy/dx)+2*(dx/dy))*T(end-jy+1,ii))+...
    k*(dy/dx)*T(end-jy+2,ii)+k*(dy/dx)*T(end-jy,ii);
    end

    % 3. Equations for bottom right-most CV
    A(ITMAX-2,1)=-k*dx/dy;
    A(ITMAX-2,2)=(k*(3*(dy/dx)+2*(dx/dy)))/alpha;
    B(ITMAX-2)=((1-alpha)/alpha)*k*(3*(dy/dx)+2*(dx/dy))*T(end-jy+1,end-1)+...
    k*(dy/dx)*T(end-jy,end-1)+k*(dy/dx)*T(end-jy+2,end-1)+2*k*(dy/dx)*T(end-jy+1,end);

    T(end-jy+1,2:end-1)=Tridiagonal(ITMAX-2,A,B);
end

% TOP-MOST
% 1. Equation for the top left-most CV
A(1,2)=(k*(3*(dy/dx)+3*(dx/dy)))/alpha;
A(1,3)=-k*dy/dx;
B(1)=((1-alpha)/alpha)*k*(3*(dy/dx)+3*(dx/dy))*T(2,2)+...
    2*k*(dy/dx)*T(2,1)+k*(dy/dx)*T(3,2)+2*k*(dy/dx)*T(1,2);

% 2. Equations for top-most internal CV's
for ii=3:JTMAX-2
    A(ii-1,1)=-k*(dx/dy);
    A(ii-1,2)=(k*(2*(dy/dx)+3*(dx/dy)))/alpha;
    A(ii-1,3)=-k*(dx/dy);
    B(ii-1)=((1-alpha)/alpha)*k*((3*(dy/dx)+2*(dx/dy))*T(2,ii))+...
    2*k*(dy/dx)*T(1,ii)+k*(dy/dx)*T(3,ii);
end

% 3. Equations for top right-most CV
A(ITMAX-2,1)=-k*dx/dy;
A(ITMAX-2,2)=(k*(3*(dy/dx)+3*(dx/dy)))/alpha;
B(ITMAX-2)=((1-alpha)/alpha)*k*(3*(dy/dx)+3*(dx/dy))*T(2,end-1)+...
    2*k*(dy/dx)*T(2,end)+k*(dy/dx)*T(3,end-1)+2*k*(dy/dx)*T(1,end-1);

T(2,2:end-1)=Tridiagonal(ITMAX-2,A,B);

% Update the right boundary
T(:,end)=-J*(dx/(2*k)) + T(:,end-1);

%% COMPUTE THE RESIDUAL
%  1. ERROR FOR EACH CV 
error=0;

%  LEFT-TOP MOST CORNER CV
error=abs((k*(3*(dy/dx)+3*(dx/dy)))*T(2,2)-(2*k*dx/dy*T(1,2)+k*dy/dx*T(2,3)+k*dx/dy*T(3,2)+2*k*dy/dx*T(2,1)));
%  LEFT-BOTTOM MOST CORNER CV
error=abs(error)+abs((k*(3*(dy/dx)+3*(dx/dy)))*T(ITMAX-1,2)-(k*dx/dy*T(ITMAX-2,2)+k*dy/dx*T(ITMAX-1,3)+2*k*dx/dy*T(ITMAX,2)+2*k*dy/dx*T(ITMAX-1,1)));
%  RIGHT-TOP MOST CORNER CV
error=abs(error)+abs((k*(3*(dy/dx)+3*(dx/dy)))*T(2,JTMAX-1)-(2*k*dx/dy*T(1,JTMAX-1)+2*k*dy/dx*T(2,JTMAX)+k*dx/dy*T(3,JTMAX-1)+k*dy/dx*T(2,JTMAX-2)));
%  RIGHT-BOTTOM MOST CORNER CV
error=abs(error)+abs((k*(3*(dy/dx)+3*(dx/dy)))*T(ITMAX-1,JTMAX-1)-(k*dx/dy*T(ITMAX-2,JTMAX-1)+2*k*dy/dx*T(ITMAX-1,JTMAX)+2*k*dx/dy*T(ITMAX,JTMAX-1)+k*dy/dx*T(ITMAX-1,JTMAX-2)));

%  2. ERROR FOR BOUNDARY CV's

%  TOP BOUNDARY
for it=3:JTMAX-2
    error=abs(error)+abs((k*(2*(dy/dx)+3*(dx/dy)))*T(2,it)-(2*k*dx/dy*T(1,it)+k*dy/dx*T(2,it+1)+k*dx/dy*T(3,it)+k*dy/dx*T(2,it-1)));
end
%  BOTTOM BOUNDARY
for ib=3:JTMAX-2
    error=abs(error)+abs((k*(2*(dy/dx)+3*(dx/dy)))*T(ITMAX-1,ib)-(k*dx/dy*T(ITMAX-2,ib)+k*dy/dx*T(ITMAX-1,ib+1)+2*k*dx/dy*T(ITMAX,ib)+k*dy/dx*T(ITMAX-1,ib-1)));
end
%  LEFT BOUNDARY
for jl=3:ITMAX-2
    error=abs(error)+abs((k*(3*(dy/dx)+2*(dx/dy)))*T(jl,2)-(k*dx/dy*T(jl-1,2)+k*dy/dx*T(jl,3)+k*dx/dy*T(jl+1,2)+2*k*dy/dx*T(jl,1)));
end
% RIGHT BOUNDARY
for jr=3:ITMAX-2
    error=abs(error)+abs((k*(3*(dy/dx)+2*(dx/dy)))*T(jr,JTMAX-1)-(k*dx/dy*T(jr-1,JTMAX-1)+2*k*dy/dx*T(jr,JTMAX)+k*dx/dy*T(jr+1,JTMAX-1)+k*dy/dx*T(jr,JTMAX-2)));
end

% 3. INTERNAL CV's
for i=3:ITMAX-2
    for j=3:JTMAX-2
        error=abs(error)+abs((k*(2*(dy/dx)+2*(dx/dy)))*T(i,j)-(k*dx/dy*T(i-1,j)+k*dy/dx*T(i,j+1)+k*dx/dy*T(i+1,j)+k*dy/dx*T(i,j-1)));
    end
end
error=abs(error);

end % end of while loop

end % end of function








