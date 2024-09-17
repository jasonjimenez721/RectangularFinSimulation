clc; clear; clf;
k = 10;   % thermal conductivity (W/m*K)
h = 50;    % convection coefficient (W/m^2*K)
Tb = 200;  % base (wall) temperature (K)
Tinf = 25; % free stream temperature (K)

dx = 0.005; % node spacing (x-direction) (m)
dy = dx;    % node spacing (y-diretion) (m)
Lx = 0.04;  % length of fin (x-dir) (m)
Ly = 0.01;  % length of fin (y-dir) (m)
Lz = 0.2 ;  %length of fin (z-dir) (m)
Nx = round(Lx/dx) + 1; % number of nodes (x-dir)
Ny = round(Ly/dy) + 1; % number of nodes (y-dir)

%For my 2D temperature matrix
T = zeros(Nx/Ny);

c1 = ((-2*h*dx)/k)*Tinf;

c2 = -2*(((h*dx)/k) + 2);  % I made shortened varibales (c's) to make the coding easier to look at and understand.

c3 = -2*(((h*dx)/k) + 1);


%   T11	T21	T31	T41	T51	T61	T71	T81	T91		T12	T22	T32	T42	T52	T62	T72	T82	T92		T13	T23	T33	T43	T53	T63	T73	T83	T93
A = [1	0	0	0	0	0	0	0	0		0	0	0	0	0	0	0	0	0		0	0	0	0	0	0	0	0	0
    0	0	0	0	0	0	0	0	0		1	0	0	0	0	0	0	0	0		0	0	0	0	0	0	0	0	0     %wall
    0	0	0	0	0	0	0	0	0		0	0	0	0	0	0	0	0	0		1	0	0	0	0	0	0	0	0
																												    
    0	1	0	0	0	0	0	0	0		1	-4	1	0	0	0	0	0	0		0	1	0	0	0	0	0	0	0
    0	0	1	0	0	0	0	0	0		0	1	-4	1	0	0	0	0	0		0	0	1	0	0	0	0	0	0
    0	0	0	1	0	0	0	0	0		0	0	1	-4	1	0	0	0	0		0	0	0	1	0	0	0	0	0
    0	0	0	0	1	0	0	0	0		0	0	0	1	-4	1	0	0	0		0	0	0	0	1	0	0	0	0    %interior
    0	0	0	0	0	1	0	0	0		0	0	0	0	1	-4	1	0	0		0	0	0	0	0	1	0	0	0
    0	0	0	0	0	0	1	0	0		0	0	0	0	0	1	-4	1	0		0	0	0	0	0	0	1	0	0
    0	0	0	0	0	0	0	1	0		0	0	0	0	0	0	1	-4	1		0	0	0	0	0	0	0	1	0
																												    
    0	0	0	0	0	0	0	0	0		0	2	0	0	0	0	0	0	0		1	c2	1	0	0	0	0	0	0
    0	0	0	0	0	0	0	0	0		0	0	2	0	0	0	0	0	0		0	1	c2	1	0	0	0	0	0
    0	0	0	0	0	0	0	0	0		0	0	0	2	0	0	0	0	0		0	0	1	c2	1	0	0	0	0
    0	0	0	0	0	0	0	0	0		0	0	0	0	2	0	0	0	0		0	0	0	1	c2	1	0	0	0   %top
    0	0	0	0	0	0	0	0	0		0	0	0	0	0	2	0	0	0		0	0	0	0	1	c2	1	0	0
    0	0	0	0	0	0	0	0	0		0	0	0	0	0	0	2	0	0		0	0	0	0	0	1	c2	1	0
    0	0	0	0	0	0	0	0	0		0	0	0	0	0	0	0	2	0		0	0	0	0	0	0	1	c2	1
																												    
    1	c2	1	0	0	0	0	0	0		0	2	0	0	0	0	0	0	0		0	0	0	0	0	0	0	0	0
    0	1	c2	1	0	0	0	0	0		0	0	2	0	0	0	0	0	0		0	0	0	0	0	0	0	0	0
    0	0	1	c2	1	0	0	0	0		0	0	0	2	0	0	0	0	0		0	0	0	0	0	0	0	0	0
    0	0	0	1	c2	1	0	0	0		0	0	0	0	2	0	0	0	0		0	0	0	0	0	0	0	0	0   %bottom
    0	0	0	0	1	c2	1	0	0		0	0	0	0	0	2	0	0	0		0	0	0	0	0	0	0	0	0
    0	0	0	0	0	1	c2	1	0		0	0	0	0	0	0	2	0	0		0	0	0	0	0	0	0	0	0
    0	0	0	0	0	0	1	c2	1		0	0	0	0	0	0	0	2	0		0	0	0	0	0	0	0	0	0
																												    
    0	0	0	0	0	0	0	0	1		0	0	0	0	0	0	0	2	c2		0	0	0	0	0	0	0	0	1    %right edge
																												    
    0	0	0	0	0	0	0	0	0		0	0	0	0	0	0	0	0	1		0	0	0	0	0	0	0	1	c3  % top corner
																												    
    0	0	0	0	0	0	0	1	c3		0	0	0	0	0	0	0	0	1		0	0	0	0	0	0	0	0	0];  % bottom corner

B = [
Tb
Tb %wall
Tb

0
0
0
0  %interior
0
0
0

c1
c1
c1
c1  %top
c1
c1
c1

c1
c1
c1
c1
c1  %bottom
c1
c1

c1 %right edge

c1  %top corner

c1]; %bottom corner

T = A\B ;


T2D(:,1) = T(1:9);
T2D(:,2) = T(10:18); % To orainze the temperatures in an orderly fashion
T2D(:,3) = T(19:27);
disp(T2D)

length = 0 : 0.005 : 0.04; %(m) %x-dir
height = 0 : 0.005: 0.01; % (m) %y-dir
[height,length] = meshgrid(height,length); %setting the plot/dimensions of the contour graph for the fin
contour(length,height,T2D, 10, 'ShowText','on'); %plots the countour line curves

Tavg = (T2D(9) + T2D(18) + T2D(27)) / 3; %Average temperature for the tip nodes
disp Tavg
disp(Tavg);

Ac = .2 * .01;
P = 2*.2 + 2*.01;
m = sqrt((h*P)/(k*Ac));
%Tb-Tinf= 175 C or 175 K because its a temperature difference
%Tinf = 25 C = 298K


Ttip = ((1/(cosh(m*Lx) + (h/(m*k)) * sinh(m*Lx))) * (175)) + 25 ;
disp Tavg1D
disp(Ttip)

M = sqrt(h*P*k*Ac) * 175;

qf1 = M * (sinh(m*Lx)+(h/(m*k))*cosh(m*Lx))/(cosh(m*Lx)+(h/(m*k))*sinh(m*Lx));
disp qf1D
disp(qf1)

q1 = (k*(dy/2*Lz)*(T(1)-T(2))/dx); %n=1
q2 = (k*(dy*Lz)*(T(10)-T(11))/dx);   %n=2
q3 = (k*(dy/2*Lz)*(T(19)-T(20))/dx); %n=3

qf2 = q1 + q2 + q3;
disp qFourier
disp(qf2)