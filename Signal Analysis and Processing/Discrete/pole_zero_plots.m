%% MATLAB script to plot the DTFT of the impulse response of an LTI system
%% knowing its transfer function H(z), through poles (2) and zeros (3)

clear all
close all
%% position of the poles 
phi=90; % degrees
p1=0.8*exp(1j*phi*pi/180.);% first pole
p1c=conj(p1);% and its complex conjugate
poles=[p1,p1c,p1,p1c];% array with all the poles
% phi2=10; % 
% p2=0.8*exp(1j*phi2*pi/180.);% another pole
% p2c=conj(p2);% and its complex conjugate
% poles=[p1,p1c,p2,p2c];% array with all the poles
%% position of the zeros (two complex conjugate plus one real, you can select all real zeros)
z1=-1.1;% first zero
zeros=[z1,1.1];% array with all the zeros
% phi3=50-30;
% z2=1*exp(1j*phi3*pi/180.);z2c=conj(z2);% another couple of zeros
% phi4=120;
% z4=1*exp(1j*phi4*pi/180.);z4c=conj(z4);% another couple of zeros
% zeros=[z1,z2,z2c,z4,z4c,1];% array with all the zeros

%% definition of the subset of the z-plane for which H(z) is evaluated
X = [-1.5:0.05:1.5];
Y = X;
[X, Y] = meshgrid(X, Y);
z=X+1j*Y;Hz=ones(size(z));
% evaluation of H(z)
for k =1: length(zeros)
    Hz = Hz.*(z-zeros(k));
end
for k =1: length(poles)
    Hz = Hz./(z-poles(k));
end
Hzabs=abs(Hz);

%% evaluation of H(e^{j2*pi*f*Delta_t})=H(e^{j*theta})
theta=[-1:0.01:1]*pi;% theta axis: 100 samples between -pi and pi
x=cos(theta); % necessary for plotting
y=sin(theta); % necessary for plotting
z=exp(1j*theta); % circumference with unit radius
Hf=ones(size(z));
for k =1: length(zeros)
    Hf = Hf.*(z-zeros(k));
end
for k =1: length(poles)
    Hf = Hf./(z-poles(k));
end
Hf=abs(Hf);
zz=Hf*0.;% zz=0 (same size as Hf) allows to plot in 3D the citrcumference of unit radius

Hfmax=max(Hf);
reflev=10;% maximum value shown in the plots is reflev*1.1
Hf1=Hf/Hfmax*reflev;
Hzabs1=min(Hzabs/Hfmax*reflev,reflev*1.10);% 

%% first 3D figure
figure()
% surface |H(z)|
mesh(X,Y,Hzabs1,...
            'facecol','no',...
            'edgealpha',.4,...
            'edgecolor','b');grid on;hold on
%contour(X,Y,Hzabs1)
% curve |H(f)| in 3D
plot3(x,y,Hf1,'g','LineWidth',2)
% circumference with radius 1 in 3D
plot3(x,y,zz,'k','LineWidth',2)
% positions of poles and zeros:
for k =1:length(poles)
    plot3(real(poles(k)),imag(poles(k)),0,'rx','LineWidth',2) % pole
end
for k =1:length(zeros)
    plot3(real(zeros(k)),imag(zeros(k)),0,'ro','LineWidth',2) % zero
end
xlabel('real(z)')
ylabel('imag(z)')
zlabel('|H(z)| (truncated peaks)')
%% second 3-D figure just to show the original peaks of H(z)
figure()
mesh(X,Y,Hzabs),hold on
plot3(x,y,zz+Hfmax,'k','LineWidth',2)
for k =1:length(poles)
    plot3(real(poles(k)),imag(poles(k)),Hfmax,'rx','LineWidth',2) % pole
end
for k =1:length(zeros)
    plot3(real(zeros(k)),imag(zeros(k)),Hfmax,'ro','LineWidth',2) % zero
end
xlabel('real(z)')
ylabel('imag(z)')
zlabel('|H(z)|')
%% plots of the DTFT of the impulse response of the system (linear y scale and dB)
figure()
plot(theta/(2*pi),Hf),grid on
xlabel('f $\Delta_t$','interpreter','latex');
ylabel('DTFT H(f)')
figure()
plot(theta/(2*pi),max(20*log10(Hf),-100)),grid on
v=axis();v(3)=max([v(3),-60]);axis(v)
xlabel('f $\Delta_t$','interpreter','latex');
ylabel('DTFT H(f) - dB')