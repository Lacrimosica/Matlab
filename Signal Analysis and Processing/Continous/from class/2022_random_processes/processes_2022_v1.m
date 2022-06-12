%% Examples of process generation
%  The random periodic process
clear all;
close all;  

% set the maximum frequency component present
T=1;    % the process period
M=10;  % maximum frequency is M/T
% the complex random variable csi_m is defiined
sigma_csi_Re=10;
sigma_csi_Im=10;
mu_csi_0=0;
sigma_csi_0=20;

% creating the process
K=50;   % oversampling factor
MK=M*K;
Dt=T/MK;   % time-step
t0=0;   % start-time for representation
t=linspace(t0,t0+(MK-1)*Dt,MK); % generating time over a single period
N_samp=MK/2+1+50;  % where the sampling instant is located in the array of MK samples
t_samp=t(N_samp);   % sampling time for the histogram


% setting the number of wanted realizations
NR=100000;
sub_sample=100;  % how many realizations to get one histogram out 
Dy=20;
y_bars=[-300:Dy:300];
Ny=length(y_bars);
pause_switch=1;
h=figure;
bars2=zeros(2*Ny,1);
y_bars2=zeros(2*Ny,1);
ticks=zeros(2,Ny);
y_ticks=zeros(2,Ny);
heights=zeros(2,Ny);
sub=0;



for nr=1:NR

% creating a single realization
csi=sigma_csi_Re*randn(M,1) + j*sigma_csi_Im*randn(M,1);   % csi_m
csi_0=randn*sigma_csi_0;
X=zeros(1,MK);
X_single=0;

if pause_switch==1
    
    % evaluating the realization over all times
    for m=1:M
        X=sqrt(1/T)*(csi(m)*exp(j*2*pi*m/T*t)+conj(csi(m))*exp(-j*2*pi*m/T*t))+X;
    end;
    X=X+csi_0;

    if nr>1
        % erases the last histogram
        set(hh,'Visible','off');
    end;

    % taking a time-sample for the histogram
    X_samp(nr)=X(N_samp);
    [bars,y_bars]=hist(X_samp,y_bars);
    % make up the histogram curve
    bars2(1:2:2*length(bars)-1)=bars;
    bars2(2:2:2*length(bars))=bars;
    y_bars2(1:2:2*length(bars)-1)=y_bars-Dy/2;
    y_bars2(2:2:2*length(bars))=y_bars+Dy/2;
    %plot3(zeros(1,MK),X,t); grid on;hold on;
    %plot(t,X); grid on;hold on;
    plot3(t,X,zeros(1,MK),'r'); grid on;hold on;
    hh=plot3(t_samp*ones(1,length(y_bars2)),y_bars2,bars2,'b','linewidth',3);
    
    % plotting the histogram ticks
    ticks(1,:)=t_samp-0.01;
    ticks(2,:)=t_samp+0.01;
    y_ticks(1,:)=y_bars+Dy/2;
    y_ticks(2,:)=y_bars+Dy/2;
    plot3(ticks,y_ticks,heights,'b','linewidth',3);
    plot3([t_samp t_samp],[y_bars(1) y_bars(end)], [ 0 0],'b','linewidth',3);
    
    title(['plotting realization number ' num2str(nr) ],'fontsize',14);
    figure(h);
    
    r=input('do you want to keep on pausing ? (y/n)   ','s');
    if r=='n'
        pause_switch=0;     
    end;
    
else
    
    % evaluating the realization only at sample time

    % evaluating the realization over all times
    for m=1:M
        X_single=sqrt(1/T)*(csi(m)*exp(j*2*pi*m/T*t_samp)+conj(csi(m))*exp(-j*2*pi*m/T*t_samp))+X_single;
    end;
    X_samp(nr)=X_single+csi_0;
    
    if mod(sub,sub_sample)==0
    
        if nr>1
           % erases the last histogram
            set(hh,'Visible','off');
        end;

        [bars,y_bars]=hist(X_samp,y_bars);
        % make up the histogram curve
        bars2(1:2:2*length(bars)-1)=bars;
        bars2(2:2:2*length(bars))=bars;
        y_bars2(1:2:2*length(bars)-1)=y_bars-Dy/2;
        y_bars2(2:2:2*length(bars))=y_bars+Dy/2;
        %plot3(zeros(1,MK),X,t); grid on;hold on;
        %plot(t,X); grid on;hold on;
        %plot3(t,X,zeros(1,MK),'r'); grid on;hold on;
        hh=plot3(t_samp*ones(1,length(y_bars2)),y_bars2,bars2,'b','linewidth',3);
        title(['plotting realization number ' num2str(nr) ]);
        figure(h);
        sub=0;
    end
    sub=sub+1;
    
end;


end;

