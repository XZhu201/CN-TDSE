function ux=difLap(ux,AB)

ks=length(ux);

Ax=AB+2;
Bx=AB-2;

Diag=zeros(1,ks);
Diag(1)=Ax;
for k=2:ks
    Diag(k)=Ax-1/Diag(k-1);
end

b=zeros(1,ks);
b(1)=Bx*ux(1)+ux(2);
b(ks)=ux(ks-1)+Bx*ux(ks);
for k=2:ks-1
    b(k)=ux(k-1)+Bx*ux(k)+ux(k+1);
end
for k=2:ks
    b(k)=b(k)+b(k-1)/Diag(k-1);
end

ux(ks)=b(ks)/Diag(ks);
for k=ks-1:-1:1
    ux(k)=(b(k)+ux(k+1))/Diag(k);
end

end