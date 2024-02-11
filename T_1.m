clear all;
close all;
clc;

%%
%naca 4 digit aerofoil (naca uvw) nomenclature
u=input('first digit of naca = ');
v=input('second digit of naca = ');
w=input('last two digit of naca = ');
ymc_c=u*0.01;
xmc_c=v*0.1;
tm=w*0.01;
U=100;
c=1;


%%
%cosine clustring
%choose even number of points on circle 
n=200;
theta0=2*pi/(n-1);
i=1:n/2;
x_c=0.5*c*(1-cos((i-0.5)*theta0));

%%
%camber line formation
for j=1:length(x_c)
    if(x_c(j)>=0) && (x_c(j)<=xmc_c)
        ycx_c(j)=ymc_c*(2*(x_c(j)/xmc_c)-(x_c(j)/xmc_c)^2);
        dyc_dx(j)=ymc_c*(2*(1/xmc_c)-(2*(x_c(j))/(xmc_c)^2));
    else
        ycx_c(j)=ymc_c*(2*((1-x_c(j))/(1-xmc_c))-((1-x_c(j))/(1-xmc_c))^2);
        dyc_dx(j)=ymc_c*(2*(-1/(1-xmc_c))+(2*(1-x_c(j))/(1-xmc_c)^2));
    end
end

%%
%
theta0=0:theta0:pi;
syms alpha
a0=alpha-(1/pi)*trapz(theta0,dyc_dx,2);
syms N
an=(2/pi)*trapz(theta0,dyc_dx.*cos(N*theta0));

b=trapz(theta0,sin(theta0).*sin(N*theta0));
h=an*b;
d=matlabFunction(h);
e=trapz(theta0,1+cos(theta0));

syms alpha
if u==0
    eqn=a0*e+d==0;
else
    eqn=a0*e+d(1)==0;
end
f=solve(eqn,alpha);
alpha=180*double(f)/pi; %in degree
zero_lift_alpha=alpha;

syms alpha x
cl=2*pi*(alpha-(1/pi)*trapz(theta0,dyc_dx.*(1-cos(theta0))));
cm_le=-cl/4+0.5*trapz(theta0,dyc_dx.*(cos(2*theta0)-cos(theta0)));
cm_qc=0.5*trapz(theta0,dyc_dx.*(cos(2*theta0)-cos(theta0)));%pitching moment at quater chord and it is independent of angle of attack
cm_x=cl*(x/c)+cm_le;
x_p=-c*cm_le/cl;


cm_le_f=matlabFunction(cm_le);
cl_f=matlabFunction(cl);
x_p_f=matlabFunction(x_p);

alpha=-pi/18:pi/180:pi/18;
cm_le_d=cm_le_f(alpha);
cl_d=cl_f(alpha);

%plots
figure(1)
yline(cm_qc);
hold on;
plot(180*alpha/pi,cl_d);
xlabel('--alpha (in degree) ---->');
ylabel('--- cm_qc and cl---->');
legend('cm_qc','cl');
hold on;
data = readtable('exp.txt');
T = data{:,1}; 
T2 = data{:,2};
plot (T,T2,'x');

figure(2)
plot(180*alpha/pi,cm_le_d);
xlabel('--alpha (in degree) ---->');
ylabel('--- cm_le ---->');

figure(3)  
plot(180*alpha/pi,x_p_f(alpha));




