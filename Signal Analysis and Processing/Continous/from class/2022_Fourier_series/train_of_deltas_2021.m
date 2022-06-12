% Train of deltas

T=1;
f=linspace(-10,10,10001);


n=0;
y=ones(1,10001);
plot(f,y);
zoom on;
grid on;
str=input('enter q to quit', 's');
h=figure;
while isempty(str)
    n=n+1;
    y=y+2*cos(2*pi*n*f*T)/T;
    plot(f,y);
    zoom on;
    grid on;
    figure(h);q
    str=input('', 's');
end;
