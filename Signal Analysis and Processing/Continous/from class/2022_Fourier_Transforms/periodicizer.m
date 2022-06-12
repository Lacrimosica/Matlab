function result=periodicizer(fun,T,n1,n2,t)

% this script takes a function fun and makes it "pseudoperiodic"


result=zeros(size(t));
for n=n1:n2
    result=fun(t-n*T)+result;
end