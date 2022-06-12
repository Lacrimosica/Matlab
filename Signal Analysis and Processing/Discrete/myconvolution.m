close all
clc
Nx=3 ;
Nh=3;
x=ones(Nx,1);
h=ones(Nh,1);
y=conv(x,h);


% % % %
my_y= zeros(Nx+Nh-1,1);
for n=0:Nx+Nh-2
    kmin=max(0,n-Nh+1);
    kmax=min(n,Nx-1);
    for k=kmin:kmax
        my_y(n+1) = my_y(n+1) +x(k+1)*h(n-k+1); 
%         we add the plus one because of indexing
    end
end

% 
figure()
n=[0:Nx+Nh-2];
stem(n,y,'o'), hold on; grid on ;
stem(n,my_y,'.');
legend('conv', 'my conv')
xlabel('n')
ylabel('y[n]')

