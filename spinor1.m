clear;
% \\\\\\\\\\\\\\\\\\\+\\\\\\\\\\\\\\\\\\\\
%            Global variables
% \\\\\\\\\\\\\\\\\\\+\\\\\\\\\\\\\\\\\\\\
global x facx kx m Nx trap c

% \\\\\\\\\\\\\\\\\\\+\\\\\\\\\\\\\\\\\\\
%              Grid setting
% \\\\\\\\\\\\\\\\\\\+\\\\\\\\\\\\\\\\\\\\

xmax = 30*pi;
Nx   = 256*8*1;
dx   = 2*xmax/Nx;
x    = -xmax+(0:Nx-1)*dx+0.5*dx;
facx = pi/xmax;
kx   = [0:Nx/2 -Nx/2+1:-1];

% \\\\\\\\\\\\\\\\\\\+\\\\\\\\\\\\\\\\\\\
%             Initail setting
% \\\\\\\\\\\\\\\\\\\+\\\\\\\\\\\\\\\\\\\\

% ==== System ====
m=1;  % mass
R=15; % width
c=10; % speed of light 

% ==== Potential ====
trapt='ho';trap=x.^2/2+disx;

pero=0;
u1=1*(pi*R^2).^(-1/4)*exp(-((x+pero).^2/(2*R^2))); %1D
u2=zeros(size(x));


if u1==u2
    u12tx='u2=u1';
elseif u2==0
    u12tx='u2=0';
end

% ==== Initial Condition ====
u0=[u1 u2];
u0=u0(:);
dt=pi/8;

% ==== Time Interval ====
tspan=0:dt:40*pi;

% ==== Simulation ====
[t,uall]=ode45(@psy1,tspan,u0);

% ==== Seperate ====
psy_p=uall(:,1:Nx);
psy_m=uall(:,Nx+1:end);

tit=['disorder-1d-spinor-spin0-' ...
    u12tx...
    '-trap=' trapt...
    '-pero=off'...
    '-xmax=' num2str(xmax) ...
    '-Nx=' num2str(Nx) ...
    '-mass=' num2str(m) ...
    '-R=' num2str(R) ...   
    '-c=' num2str(c) ...
    '-disdump=' num2str(A)...
    '-tspan=' num2str(round(tspan(1)/pi)) '-' num2str(round(tspan(end)/pi)) ];

tit1=['save ' tit '.mat'];

% ==== Save ====
mat_t=[tit '.mat'];
eval(tit1);

% \\\\\\\\\\\\\\\\\\\+\\\\\\\\\\\\\\\\\\\
%                Equation
% \\\\\\\\\\\\\\\\\\\+\\\\\\\\\\\\\\\\\\\\
function dphidt=psy1(t,u)
global x kx facx m Nx trap c

% ==== Determinate ====
psy_p=u(1:Nx); psy_p=reshape(psy_p,size(x));
psy_m=u(Nx+1:end); psy_m=reshape(psy_m,size(x));

mx=fft(abs(psy_p).^2-abs(psy_m).^2)*(x(2)-x(1));
mx=mx(1);

tr=sprintf('t= %0.4e mx= %0.4e',t,mx);
disp(tr);

% ==== Eq ====
ug=psy_m+psy_p;
uh=fft(ug);

dudx2=ifft(-uh.*kx.^2)*facx^2;

dpsydt_p = -1i*(-dudx2/(2)+c.^2*m*psy_p+trap.*psy_p);
dpsydt_m = -1i*( dudx2/(2)-c.^2*m*psy_m+trap.*psy_m);

dphidt=[dpsydt_p(:) dpsydt_m(:)];

dphidt=dphidt(:);

end