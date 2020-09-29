n=100;
s=ones(1,n-1);
A=eye(n)*10+diag(s,-1)+diag(s,1);
ty=ones(1,n-2)*12;
b=rand(100,1)*5;
t0=cputime;
for k=1:1:n
    A(k,k)=sqrt(A(k,k))
end
for k=1:1:n-1
    A(k+1:n,k)=A(k+1:n,k)/A(k,k)
    for j=k+1:1:n
        A(j:n,j)=A(j:n,j)-A(j:n,k)*A(j,k)
    end
end
L=tril(A)
y=zeros(n,1)
y(1)=b(1)/L(1,1)
for k=2:1:n
    y(k)=(b(k)-dot(L(k,1:k-1),y(1:k-1)))/L(k,k)
end
x=zeros(n,1)
x(n)=y(n)/L(n,n)
for k=n-1:-1:1
    x(k)=(y(k)-dot(L(k+1:n,k),x(k+1:n)))/L(k,k);
end
t=cputime-t0


