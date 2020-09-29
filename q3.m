A=[0.05,0.07,0.06,0.05;0.07,0.10,0.08,0.07;0.06,0.08,0.10,0.09;0.05,0.07,0.09,0.10];
b=[0.23;0.32;0.33;0.31];
n=4
A1=A
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
    x(k)=(y(k)-L(k+1:n)*x(k+1:n))/L(k,k);
end
