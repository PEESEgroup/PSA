function [results, parameters]=sortt(A)
%% 

a=A.pops(end, :);

[~, n]=size(a);

q=length(a(1, 1).obj);
results=zeros(n, q);
qq=a(1, 1).var;
[~, nn]=size(qq);

parameters=zeros(n, nn);


for i= 1:n
    
   results(i, 1:q)=a(1, i).obj;
   parameters(i, 1:nn)=a(1, i).var;
   
end
%% 

results=abs(results);

end