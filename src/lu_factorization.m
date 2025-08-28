function [U,L, U_int, L_int]= lu_factorization(A,A_int)

[m,n]=size(A);
U=zeros(m);
L=zeros(m);

U_int=int16(zeros(m));
L_int=int16(zeros(m));

for j=1:m
	U(1,j)=A(1,j);
    U_int(1,j)=A_int(1,j);
end

for j=1:m
	L(j,j)=1;
    L_int(j,j)=1;
end

for i=2:m
    %% L
    for k=1:i-1
        prod_L=0;
        prod_L_int=int16(0);
        if k==1
            prod_L=0;
            prod_L_int=0;
        else
            for p=1:k-1
                prod_L=prod_L+L(i,p)*U(p,k);
                prod_L_int=int32(prod_L_int)+int32(L_int(i,p))*int32(U_int(p,k));
                prod_L_int = int16(round(prod_L_int/65536));
            end
        end
        L(i,k)=(A(i,k)-prod_L)/U(k,k);
        L_int(i,k) = cordic_divide(A_int(i,k)-prod_L_int, U_int(k,k), 11);
    end
    %% U
    for k=i:m
        prod_U=0;
        prod_U_int=int16(0);
        for p=1:i-1
            prod_U=prod_U+L(i,p)*U(p,k);
            prod_U_int=int32(prod_U_int)+int32(L_int(i,p))*int32(U_int(p,k));
            prod_U_int = int16(round(prod_U_int/65536));
        end
        U(i,k)=A(i,k)-prod_U;
        U_int(i,k)=A_int(i,k)-prod_U_int;
        end
    end
end