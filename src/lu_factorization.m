function [U,L, U_int, L_int]= lu_factorization(A,A_int)

A_fi = fi(A_int,1,12,0);

[m,n]=size(A);
U=zeros(m);
L=zeros(m);

U_fi=fi(zeros(m),1,12,0);
L_fi=fi(zeros(m),1,15,0);

for j=1:m
	U(1,j)=A(1,j);
    U_fi(1,j)=A_fi(1,j);
end

for j=1:m
	L(j,j)=1;
    L_fi(j,j)=1;
end

for i=2:m
    %% L
    for k=1:i-1
        prod_L=0;
        prod_L_fi=fi(0,1,12,0);
        if k==1
            prod_L=0;
            prod_L_fi=fi(0,1,12,0);
        else
            for p=1:k-1

                b = L(i,p)*U(p,k);
                prod_L = prod_L + b;

                b1 = (L_fi(i,p)) * (U_fi(p,k)); % fi(1,15,11) + fi(1,12,5) = fi(1,27,16)
                b1 = fi(bitshift(int32(b1),-6,'int32'),1,21,0); % fi(1,21,10)
                prod_L_fi = fi(prod_L_fi + b1,1,18,0); % fi(1,12,10) + fi(1,21,10) = fi(1,22,10)
                % prod_L_fi = fi(bitshift(int32(prod_L_fi),-6,'int32'),1,16,0); % fi(1,16,4)

            end
        end

        if i == 3 && k == 2
            q = 1;
        end

        dd = A(i,k) - prod_L;
        rr = dd / U(k,k);
        L(i,k) = rr;

        num = fi(A_fi(i,k) - prod_L_fi,1,13,0); % fi(1,12,10) - fi(1,12,10) = fi(1,13,10)

        % num_bit = dec2bin(num);
        % l = length(num_bit);
        % 
        % if (num_bit(l) == true)
        %     num = fi(bitshift(int16(num),-1,'int16'),1,12,0); % fi(1,12,10) % round
        %     num = num + 1;
        % else
        %     num = fi(bitshift(int16(num),-1,'int16'),1,12,0); % fi(1,12,9)
        % end

        double(num)*2^-10

        aa = fi(cordic_divide(num, U_fi(k,k), 11),1,16,0);  % fi(1,16,11)
        L_fi(i,k) = aa; % fi(1,16,11)
        L_fi_d(i,k) = aa * 2^-11;
    end
    %% U
    for k=i:m
        prod_U=0;
        prod_U_fi=fi(0,1,12,0);
        for p=1:i-1
            a = L(i,p)*U(p,k);
            prod_U = prod_U + a;
            %%

            a1 = (L_fi(i,p)) * (U_fi(p,k)); % fi(1,16,11) * fi(1,12,10) = fi(1,28,21)
            a1 = fi(bitshift(int32(a1),-11,'int32'),1,16,0); % fi(1,17,10)

            prod_U_fi = fi(prod_U_fi + a1,1,17,0); % fi(1,12,10) + fi(1,17,10) = fi(1,18,10)
            % prod_U_fi = fi(bitshift(int32(prod_U_fi),-2,'int32'),1,16,0); % fi(1,16,8)

        end
        U(i,k) = A(i,k) - prod_U;


        if (i == 2 && k == 2)
            e = 1;
        end

        ty = A_fi(i,k) - prod_U_fi; % fi(1,12,10) - (1,18,10) = 1,19,10
        ty = fi(bitshift(int32(ty),-7,'int32'),1,12,0); % fi(1,12,3)

        U_fi(i,k) = ty; 
        end
    end
end