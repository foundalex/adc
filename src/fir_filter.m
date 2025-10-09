function y = fir_filter(b, x, width_in, width_out)

    z = fi(zeros(size(b)),1,width_in,0);
    y = fi(zeros(size(x)),1,width_out,0);

    p = 0;
    nx = length(x);
    nb = length(b);



    for n=1:nx
        p=p+1; if p>nb, p=1; end
        z(p) = x(n);
        acc = fi(0,1,width_out,0);
        k = p;
        for j=1:nb
            acc(:) = acc + (b(j,:)) * (z(k));
            k=k-1; if k<1, k=nb; end
        end
        y(n) = acc;
    end

end
