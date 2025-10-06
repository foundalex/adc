function y = fir_filter(b,x)

    z = int32(zeros(size(b)));
    y = int32(zeros(size(x)));

    p = 0;
    nx = length(x);
    nb = length(b);



    for n=1:nx
        p=p+1; if p>nb, p=1; end
        z(p) = x(n);
        acc = int32(0);
        k = p;
        for j=1:nb
            acc = acc + int32(b(j,:)) * (z(k));
            k=k-1; if k<1, k=nb; end
        end
        y(n) = acc;
    end

end
