N=10
A1=zeros(10,10)
for j=1:1:10
    for i=j:1:10
    A1(i,j)=(i+j)^2
    end
end
B=eye(10)
A=[A1 B]

 for k=1:1:N
   A(k,:)=A(k,:)/A(k,k)
   A(k+1:10,k:20)=A(k+1:10,k:20)-A(k+1:10,k)*A(k,k:20)
 end
B=A(1:10,11:20)