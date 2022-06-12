
%%  Examples of pdf estimation of single and joint RVs

%% single uniformly distributed RV
close all;
clear all;

N=100000;  % how many experiments
csi=rand(1,N);   % the values of csi

% plotting the histogram
Dx=0.1;    % box width (i.e., value of Delta_x)
[box,x]=hist(csi,[-1+Dx/2:Dx:2-Dx/2]);
figure;bar(x,box/N,1);grid on; 
title(['raw histogram, Dx = ' num2str(Dx) ' , boxes sum to ' num2str(sum(box)/N) ]);
pause;

% plotting the pdf estimate
pdf=box/N/Dx;
figure;bar(x,pdf,1);grid on;
title(['pdf estimate (integrates to ' num2str(sum(pdf)*Dx) ')']);
pause;

% plotting the expected analytical pdf
ax=[min(x):0.001:max(x)];
apdf=HPi(1,ax-1/2);
hold on; plot(ax,apdf,'r--', 'linewidth',2);

return;

%% single Gaussian-distributed RV
clear all;

N=10000;  % how many experiments
sigma=1;
mu=0;
csi=sigma*randn(1,N)+mu;   % the values

% plotting the histogram
Dx=0.1*sigma;    % box width (i.e., value of Delta_x)
x=[-6*sigma:Dx:6*sigma];
[box,x]=hist(csi,x);
figure;bar(x,box/N,1);grid on; 
title(['raw histogram, Dx = ' num2str(Dx) ' , boxes sum to ' num2str(sum(box)/N) ]);
pause;

% plotting the pdf estimate
pdf=box/N/Dx;
figure;bar(x,pdf,1);grid on;
title(['pdf estimate (integrates to ' num2str(sum(pdf)*Dx) ')']);
pause;

% plotting the expected analytical pdf
apdf=1/sqrt(2*pi*sigma^2)*exp(-(x-mu).^2 /2/(sigma^2));
hold on; plot(x,apdf,'r--', 'linewidth',2);

return;

%% Rayleigh-distributed RV
clear all;

N=1000000;  % how many experiments
% variance of the component Gaussian variables
sigma=1;
%mu=2;
csi=sqrt((sigma*randn(1,N)).^2+(sigma*randn(1,N)).^2);   % the values

% mean value and standard deviation of the Rayleigh
mean_value=sqrt(pi/2)*sigma;
std_dev=sqrt((4-pi))/2*sigma;
% plotting the histogram
Dx=0.01*sigma;    % box width (i.e., value of Delta_x)

x=[-2*std_dev:Dx:mean_value+8*std_dev];
[box,x]=hist(csi,x);
figure;bar(x,box/N,1);grid on; 
title(['raw histogram, Dx = ' num2str(Dx) ' , boxes sum to ' num2str(sum(box)/N) ]);
pause;

% plotting the pdf estimate
pdf=box/N/Dx;
figure;bar(x,pdf,1);grid on;
title(['pdf estimate (integrates to ' num2str(sum(pdf)*Dx) ')']);
pause;

% plotting the expected analytical pdf
apdf=x./sigma^2.*exp(-x.^2/(2*sigma^2)).*(x>0);
%1/sqrt(2*pi*sigma^2)*exp(-(x-mu).^2 /2/(sigma^2));
hold on; plot(x,apdf,'r--', 'linewidth',2);

return;

%% Independent jointly-distributed Gaussian RVs
clear all;

N=1000000;  % how many experiments
sigma_csi=10;
sigma_eta=10;
mu_csi=15;
mu_eta=15;
csi_eta=zeros(N,2);
csi_eta(:,1)=sigma_csi*randn(1,N)+mu_csi;   % csi: temperature in Torino
csi_eta(:,2)=sigma_eta*randn(1,N)+mu_eta;   % eta: temperature in Beijing

% plotting the histogram
Dx=0.25*sigma_csi;    % box width (i.e., value of Delta_x)
Dy=0.25*sigma_eta;    % box width (i.e., value of Delta_y)

reticulate={ (-6*sigma_csi+mu_csi):Dx:(6*sigma_csi+mu_csi)  (-6*sigma_eta+mu_eta):Dy:(6*sigma_eta+mu_eta) };
figure;hist3(csi_eta, reticulate);[box]=hist3(csi_eta, reticulate);
title(['raw histogram: number of occurrences per box, Dx=' num2str(Dx) ', Dy=' num2str(Dx)]);
pause;

csi_eta_pdf=box/Dx/Dy/N;    % estimating the pdf
figure; mesh(reticulate{1},reticulate{2},csi_eta_pdf);axis square;
title(['estimate of joint pdf of \xi and \eta - integrates to ' num2str(sum(sum(csi_eta_pdf))*Dx*Dy)]);
pause;

figure;contour(reticulate{1},reticulate{2},csi_eta_pdf);axis square;grid on;
title(['estimate of joint pdf of \xi and \eta - integrates to ' num2str(sum(sum(csi_eta_pdf))*Dx*Dy)]);
pause;



%% Transforming a uniform RV
clear all;

N=10000000;  % how many experiments
extension=pi;
mu=0;
phi=extension*(rand(1,N)-1/2-mu);   % zero-mean with given extension

% plotting the pdf estimate
Dx=0.02*extension;    % box width (i.e., value of Delta_x)
[box,x]=hist(phi,[-extension+Dx/2:Dx:extension-Dx/2]);
pdf=box/N/Dx;
figure;bar(x,pdf);grid on;
title(['estimate of pdf of phi - (integrates to ' num2str(sum(pdf)*Dx) ')']);

pause;

% transforming the variable
eta=cos(phi);

% plotting the pdf estimate
Dy=0.05;    % box width (i.e., value of Delta_x)
[box,y]=hist(eta,[-1+Dy/2:Dy:2-Dy/2]);
pdf_eta=box/N/Dy;
figure;bar(y,pdf_eta);grid on;
title(['estimate of pdf of \eta=cos(\phi) - (integrates to ' num2str(sum(pdf_eta)*Dy) ')']);
return;


%% sum of independent uniform RVs and central limit theorem

clear all;

N=1000000;  % how many experiments
extension=sqrt(12);    % extension of uniform pdf
mu=0;   % mean value of uniform pdf

eta=zeros(1,N);     % initializing array

for k=1:100

csi=extension*(rand(1,N)-1/2-mu);   % creating a new uniformly-distributed RV
eta=eta+csi;    % summing to the previous result

% plotting the pdf estimate
Dx=0.01*extension*sqrt(k);    % box width (i.e., value of Delta_x)
[box,x]=hist(eta,[-1.5*sqrt(k)*extension:Dx:sqrt(k)*extension*1.5]);
pdf=box/N/Dx;
figure;bar(x,pdf,1);grid on;

var_eta=var(eta);   % estimating variance of eta
mean_eta=mean(eta); % estimating mean of eta
title(['estimate of pdf of \eta - integrates to ' num2str(sum(pdf)*Dx) ' - mean \eta=' num2str(mean_eta) ' - var \eta=' num2str(var_eta)]);

% plotting the expected analytical pdf
sigma_eta=extension/2/sqrt(3)*sqrt(k);
mu_eta=k*mu;
apdf=1/sqrt(2*pi*sigma_eta^2)*exp(-(x-mu_eta).^2 /2/(sigma_eta^2));
hold on; plot(x,apdf,'r--', 'linewidth',2);hold off;

pause;

end % index k

pause;


%% Independent jointly distributed uniform RVs

clear all;

N=100000;  % how many experiments
sigma_csi=3;
sigma_eta=3;
mu_csi=2;
mu_eta=4;
csi_eta=zeros(N,2); % first column contains csi, second column contains eta
csi_eta(:,1)=sigma_csi*sqrt(12)*(rand(1,N)-1/2)+mu_csi;   % the values of csi
csi_eta(:,2)=sigma_eta*sqrt(12)*(rand(1,N)-1/2)+mu_eta;   % the values of eta

% plotting the histogram
Dx=0.2*sigma_csi;    % box width (i.e., value of Delta_x)
Dy=0.2*sigma_eta;    % box width (i.e., value of Delta_y)

reticulate={ (-6*sigma_csi+mu_csi):Dx:(6*sigma_csi+mu_csi)  (-6*sigma_eta+mu_eta):Dy:(6*sigma_eta+mu_eta) };
figure;hist3(csi_eta, reticulate);[box]=hist3(csi_eta, reticulate);
title(['raw histogram: number of occurrences per box, Dx=' num2str(Dx) ', Dy=' num2str(Dx)]);
pause;

csi_eta_pdf=box/Dx/Dy/N;    % estimating the pdf
figure; mesh(reticulate{1},reticulate{2},csi_eta_pdf);
title(['estimate of joint pdf of \xi and \eta - integrates to ' num2str(sum(sum(csi_eta_pdf))*Dx*Dy)]);
pause;

figure;contour(reticulate{1},reticulate{2},csi_eta_pdf);axis square;grid on;
title(['estimate of joint pdf of \xi and \eta - integrates to ' num2str(sum(sum(csi_eta_pdf))*Dx*Dy)]);
pause;


%% Dependent jointly distributed uniform RVs

clear all;

N=1000000;  % how many experiments
w_csi=1;            % FWHM of the distribution of csi
w_rho=1;            % FWHM of the distribution of rho
sigma_csi=w_csi/2/sqrt(3);
sigma_rho=w_rho/2/sqrt(3);
mu_csi=0;
mu_rho=0;
csi_eta=zeros(N,2); % first column contains csi, second column contains eta
csi=sigma_csi*sqrt(12)*(rand(1,N)-1/2)+mu_csi;   % the values of csi
rho=sigma_rho*sqrt(12)*(rand(1,N)-1/2)+mu_rho;   % the values of eta
csi_eta(:,1)=csi;
csi_eta(:,2)=csi+rho;
% plotting the histogram
Dx=0.1*sigma_csi;    % box width (i.e., value of Delta_x)
Dy=0.1*sigma_csi;    % box width (i.e., value of Delta_y)

reticulate={ (-6*sigma_csi+mu_csi):Dx:(6*sigma_csi+mu_csi)  (-6*sigma_csi+mu_csi):Dy:(6*sigma_csi+mu_csi) };
figure;hist3(csi_eta, reticulate);[box]=hist3(csi_eta, reticulate);
title(['raw histogram: number of occurrences per box, Dx=' num2str(Dx) ', Dy=' num2str(Dx)]);
pause;

csi_eta_pdf=box/Dx/Dy/N;    % estimating the pdf
figure; mesh(reticulate{1},reticulate{2},csi_eta_pdf);
title(['estimate of joint pdf of \xi and \eta - integrates to ' num2str(sum(sum(csi_eta_pdf))*Dx*Dy)]);
pause;

figure;contour(reticulate{1},reticulate{2},csi_eta_pdf);axis square;grid on;
title(['estimate of joint pdf of \xi and \eta - integrates to ' num2str(sum(sum(csi_eta_pdf))*Dx*Dy)]);
pause;


%% Independent jointly-distributed Gaussian RVs
clear all;

N=1000000;  % how many experiments
sigma_csi=10;
sigma_eta=10;
mu_csi=15;
mu_eta=15;
csi_eta=zeros(N,2);
csi_eta(:,1)=sigma_csi*randn(1,N)+mu_csi;   % csi: temperature in Torino
csi_eta(:,2)=sigma_eta*randn(1,N)+mu_eta;   % eta: temperature in Beijing

% plotting the histogram
Dx=0.25*sigma_csi;    % box width (i.e., value of Delta_x)
Dy=0.25*sigma_eta;    % box width (i.e., value of Delta_y)

reticulate={ (-6*sigma_csi+mu_csi):Dx:(6*sigma_csi+mu_csi)  (-6*sigma_eta+mu_eta):Dy:(6*sigma_eta+mu_eta) };
figure;hist3(csi_eta, reticulate);[box]=hist3(csi_eta, reticulate);
title(['raw histogram: number of occurrences per box, Dx=' num2str(Dx) ', Dy=' num2str(Dx)]);
pause;

csi_eta_pdf=box/Dx/Dy/N;    % estimating the pdf
figure; mesh(reticulate{1},reticulate{2},csi_eta_pdf);axis square;
title(['estimate of joint pdf of \xi and \eta - integrates to ' num2str(sum(sum(csi_eta_pdf))*Dx*Dy)]);
pause;

figure;contour(reticulate{1},reticulate{2},csi_eta_pdf);axis square;grid on;
title(['estimate of joint pdf of \xi and \eta - integrates to ' num2str(sum(sum(csi_eta_pdf))*Dx*Dy)]);
pause;


%% Partially dependent jointly-distributed Gaussian RVs

clear all;

N=1000000;  % how many experiments
sigma_rho=10;
sigma_delta1=3;
sigma_delta2=3;
mu_rho=15;
mu_delta1=0;
mu_delta2=0;
rho=zeros(N,1);
delta1=zeros(N,1);
delta2=zeros(N,1);
csi_eta=zeros(N,2);

rho=sigma_rho*randn(1,N)+mu_rho;   % the values
delta1=sigma_delta1*randn(1,N)+mu_delta1;   % the values
delta2=sigma_delta2*randn(1,N)+mu_delta2;   % the values

csi_eta(:,1)=rho+delta1;    % csi, temperature in Torino
csi_eta(:,2)=rho+delta2;    % eta, temperature in Moncalieri

% plotting the histogram
Dx=0.2*sigma_rho;    % box width (i.e., value of Delta_x)
Dy=0.2*sigma_rho;    % box width (i.e., value of Delta_y)

reticulate={ (-6*sigma_rho+mu_rho):Dx:(6*sigma_rho+mu_rho)  (-6*sigma_rho+mu_rho):Dy:(6*sigma_rho+mu_rho) };
figure;hist3(csi_eta, reticulate);[box]=hist3(csi_eta, reticulate);
title(['raw histogram: number of occurrences per box, Dx=' num2str(Dx) ', Dy=' num2str(Dx)]);
pause;

csi_eta_pdf=box/Dx/Dy/N;    % estimating the pdf
cc=corrcoef(csi_eta(:,1),csi_eta(:,2)); % correlation coefficient
cov(csi_eta(:,1),csi_eta(:,2)); % covariance matrix
figure; mesh(reticulate{1},reticulate{2},csi_eta_pdf);axis square;
title(['estimate of joint pdf of \xi and \eta - integrates to ' num2str(sum(sum(csi_eta_pdf))*Dx*Dy) ' - corr. coeff.: ' num2str(cc(1,2))]);
pause;

figure;contour(reticulate{1},reticulate{2},csi_eta_pdf);axis square;grid on;
title(['estimate of joint pdf of \xi and \eta - integrates to ' num2str(sum(sum(csi_eta_pdf))*Dx*Dy) ' - corr. coeff.: ' num2str(cc(1,2))]);
pause;

% %% Fully dependent distributed Gaussian RVs
% 
% csi_eta(:,2)=csi_eta(:,1);
% figure;hist3(csi_eta, reticulate);[box]=hist3(csi_eta, reticulate);
% pause;
% csi_eta_pdf=box/Dx/Dy/N; 
% figure; mesh(reticulate{1},reticulate{2},csi_eta_pdf);
% pause;
% figure;contour(reticulate{1},reticulate{2},csi_eta_pdf);axis square;grid on;

%% Partially dependent jointly-distributed Gaussian RVs *analytical* (no histogram)

clear all;

mu_csi=0;
mu_eta=0;

sigma_csi=1;
sigma_eta=1;

rho_csi_eta=0.7;


% plotting the pdf
Dx=0.2*sigma_csi;    % box width (i.e., value of Delta_x)
Dy=0.2*sigma_eta;    % box width (i.e., value of Delta_y)

[x y]=meshgrid( (-6*sigma_csi+mu_csi):Dx:(6*sigma_csi+mu_csi) , (-6*sigma_eta+mu_eta):Dy:(6*sigma_eta+mu_eta) );
%
figure; mesh(x,y,1/sqrt(2*pi*sigma_csi^2)*1/sqrt(2*pi*sigma_eta^2)*1/sqrt(1-rho_csi_eta^2)...
    .* exp( ( - (x-mu_csi).^2 ./ 2 ./ sigma_csi.^2  - (y-mu_eta).^2 ./ 2 ./ sigma_eta.^2 ...
    +rho_csi_eta .* (x-mu_csi) ./ sigma_csi .* (y-mu_eta) ./ sigma_eta )...
    ./(1-rho_csi_eta^2) )... 
    );
%     .* exp( ( - (x-mu_csi).^2 ./ 2 ./ sigma_csi.^2  - (y-mu_eta).^2 ./ 2 ./ sigma_eta.^2 +...
%     rho_csi_eta .* (x-mu_csi) ./ sigma_csi .* (y-mu_eta) ./ sigma_eta  ) ./(1-rho_csi_eta^2)  );
axis square;
title(['joint pdf of \xi and \eta - corr. coeff.: ' rho_csi_eta]);
h=gca;
h.XAxis.FontSize=16;
h.YAxis.FontSize=16;
h.ZAxis.FontSize=16;
pause;

figure;contour(x,y,...
    1/sqrt(2*pi*sigma_csi^2)*1/sqrt(2*pi*sigma_eta^2)*1/sqrt(1-rho_csi_eta^2)...
    .* exp( ( - (x-mu_csi).^2 ./ 2 ./ sigma_csi.^2  - (y-mu_eta).^2 ./ 2 ./ sigma_eta.^2 ...
    +rho_csi_eta .* (x-mu_csi) ./ sigma_csi .* (y-mu_eta) ./ sigma_eta )...
    ./(1-rho_csi_eta^2) )... 
);axis square;grid on;
title(['joint pdf of \xi and \eta - corr. coeff.: ' rho_csi_eta]);
h=gca;
h.XAxis.FontSize=16;
h.YAxis.FontSize=16;

pause;


% %% Fully dependent distributed Gaussian RVs
% 
% csi_eta(:,2)=csi_eta(:,1);
% figure;hist3(csi_eta, reticulate);[box]=hist3(csi_eta, reticulate);
% pause;
% csi_eta_pdf=box/Dx/Dy/N; 
% figure; mesh(reticulate{1},reticulate{2},csi_eta_pdf);
% pause;
% figure;contour(reticulate{1},reticulate{2},csi_eta_pdf);axis square;grid on;


%% Examples of process generation
%  The random periodic process

% set the maximum frequency component present
T=1;    % the process period
M=100;  % maximum frequency is M/T
% the complex random variable csi_m is defiined
sigma_csi_Re=10;
sigma_csi_Im=10;
mu_csi_0=0;
sigma_csi_0=10;

% creating the process
K=100;   % oversampling factor
MK=M*K;
Dt=T/MK;   % time-step
t0=0;   % start-time for representation
t=linspace(t0,t0+(MK-1)*Dt,MK); % generating time over a single period

% setting the number of wanted realizations
NR=1;

figure;

for nr=1:NR

% creating a single realization
X=zeros(1,MK);
csi=sigma_csi_Re*randn(1,M) + j*sigma_csi_Im*randn(1,M);   % csi_m
csi_0=randn*sigma_csi_0;

for m=1:M
    X=sqrt(1/T)*(csi(m)*exp(j*2*pi*m/T*t)+conj(csi(m))*exp(-j*2*pi*m/T*t))+X;
end;
    X=X+csi_0;

plot(t,X); grid on;hold on;

end;

return





