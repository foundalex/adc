function [y, y1_shift_round] = fir_filter(b, x, width_in, width_out)

    z = fi(zeros(size(b)),1,width_in,0);
    y = fi(zeros(size(x)),1,width_out,0);

    p = 0;
    nx = length(x);
    nb = length(b);

    buffer = int32(zeros(1,length(b)));

    % shift = 0;
    % shift_sum = 0;

    coeff_int = int32(b);

	mult2 = int32(zeros(length(x),1));
    mult3 = int32(zeros(length(x),1));
    mult4 = int32(zeros(length(x),1));
    mult5 = int32(zeros(length(x),1));
	mult6 = int32(zeros(length(x),1));
	mult8 = int32(zeros(length(x),1));
	mult7 = int32(zeros(length(x),1));
	mult9 = int32(zeros(length(x),1));
	mult10 = int32(zeros(length(x),1));
	mult11 = int32(zeros(length(x),1));
	mult12 = int32(zeros(length(x),1));
	mult13 = int32(zeros(length(x),1));
	mult14 = int32(zeros(length(x),1));
	mult15 = int32(zeros(length(x),1));
	mult16 = int32(zeros(length(x),1));
	mult17 = int32(zeros(length(x),1));
	mult18 = int32(zeros(length(x),1));
	mult19 = int32(zeros(length(x),1));
	mult20 = int32(zeros(length(x),1));
	mult21 = int32(zeros(length(x),1));
	mult22 = int32(zeros(length(x),1));
	mult23 = int32(zeros(length(x),1));
	mult24 = int32(zeros(length(x),1));
	mult25 = int32(zeros(length(x),1));
	mult26 = int32(zeros(length(x),1));
	mult27 = int32(zeros(length(x),1));
	mult28 = int32(zeros(length(x),1));
	mult29 = int32(zeros(length(x),1));
	mult30 = int32(zeros(length(x),1));
	mult31 = int32(zeros(length(x),1));
	mult32 = int32(zeros(length(x),1));
	mult33 = int32(zeros(length(x),1));
	mult34 = int32(zeros(length(x),1));
	mult35 = int32(zeros(length(x),1));
	mult36 = int32(zeros(length(x),1));
	mult37 = int32(zeros(length(x),1));
	mult38 = int32(zeros(length(x),1));
	mult39 = int32(zeros(length(x),1));
	mult40 = int32(zeros(length(x),1));
	mult41 = int32(zeros(length(x),1));
	mult42 = int32(zeros(length(x),1));
	mult43 = int32(zeros(length(x),1));
	mult44 = int32(zeros(length(x),1));
	mult45 = int32(zeros(length(x),1));
	mult46 = int32(zeros(length(x),1));
	mult47 = int32(zeros(length(x),1));
	mult48 = int32(zeros(length(x),1));
	mult49 = int32(zeros(length(x),1));
	mult50 = int32(zeros(length(x),1));
	mult51 = int32(zeros(length(x),1));
	mult52 = int32(zeros(length(x),1));
	mult53 = int32(zeros(length(x),1));
	mult54 = int32(zeros(length(x),1));
	mult55 = int32(zeros(length(x),1));
	mult56 = int32(zeros(length(x),1));
	mult57 = int32(zeros(length(x),1));
	mult58 = int32(zeros(length(x),1));
	mult59 = int32(zeros(length(x),1));
	mult60 = int32(zeros(length(x),1));
	mult61 = int32(zeros(length(x),1));
	mult62 = int32(zeros(length(x),1));
	mult63 = int32(zeros(length(x),1));
	mult64 = int32(zeros(length(x),1));
	mult65 = int32(zeros(length(x),1));
	mult66 = int32(zeros(length(x),1));
	mult67 = int32(zeros(length(x),1));
	mult68 = int32(zeros(length(x),1));
	mult69 = int32(zeros(length(x),1));
	mult70 = int32(zeros(length(x),1));
	mult71 = int32(zeros(length(x),1));
	mult72 = int32(zeros(length(x),1));
	
    mult_overflow2 = (zeros(length(x),1));
	mult_overflow3 = (zeros(length(x),1));
	mult_overflow4 = (zeros(length(x),1));
	mult_overflow5 = (zeros(length(x),1));
	mult_overflow6 = (zeros(length(x),1));
	mult_overflow7 = (zeros(length(x),1));
	mult_overflow8 = (zeros(length(x),1));
	mult_overflow9 = (zeros(length(x),1));
	mult_overflow10 = (zeros(length(x),1));
	mult_overflow11 = (zeros(length(x),1));
	mult_overflow12 = (zeros(length(x),1));
	mult_overflow13 = (zeros(length(x),1));
	mult_overflow14 = (zeros(length(x),1));
	mult_overflow15 = (zeros(length(x),1));
	mult_overflow16 = (zeros(length(x),1));
	mult_overflow17 = (zeros(length(x),1));
	mult_overflow18 = (zeros(length(x),1));
	mult_overflow19 = (zeros(length(x),1));
	mult_overflow20 = (zeros(length(x),1));
	mult_overflow21 = (zeros(length(x),1));
	mult_overflow22 = (zeros(length(x),1));
	mult_overflow23 = (zeros(length(x),1));
	mult_overflow24 = (zeros(length(x),1));
	mult_overflow25 = (zeros(length(x),1));
	mult_overflow26 = (zeros(length(x),1));
	mult_overflow27 = (zeros(length(x),1));
	mult_overflow28 = (zeros(length(x),1));
	mult_overflow29 = (zeros(length(x),1));
	mult_overflow30 = (zeros(length(x),1));
	mult_overflow31 = (zeros(length(x),1));
	mult_overflow32 = (zeros(length(x),1));
	mult_overflow33 = (zeros(length(x),1));
	mult_overflow34 = (zeros(length(x),1));
	mult_overflow35 = (zeros(length(x),1));
	mult_overflow36 = (zeros(length(x),1));
	mult_overflow37 = (zeros(length(x),1));
	mult_overflow38 = (zeros(length(x),1));
	mult_overflow39 = (zeros(length(x),1));
	mult_overflow40 = (zeros(length(x),1));
	mult_overflow41 = (zeros(length(x),1));
	mult_overflow42 = (zeros(length(x),1));
	mult_overflow43 = (zeros(length(x),1));
	mult_overflow44 = (zeros(length(x),1));
	mult_overflow45 = (zeros(length(x),1));
	mult_overflow46 = (zeros(length(x),1));
	mult_overflow47 = (zeros(length(x),1));
	mult_overflow48 = (zeros(length(x),1));
	mult_overflow49 = (zeros(length(x),1));
	mult_overflow50 = (zeros(length(x),1));
	mult_overflow51 = (zeros(length(x),1));
	mult_overflow52 = (zeros(length(x),1));
	mult_overflow53 = (zeros(length(x),1));
	mult_overflow54 = (zeros(length(x),1));
	mult_overflow55 = (zeros(length(x),1));
	mult_overflow56 = (zeros(length(x),1));
	mult_overflow57 = (zeros(length(x),1));
	mult_overflow58 = (zeros(length(x),1));
	mult_overflow59 = (zeros(length(x),1));
	mult_overflow60 = (zeros(length(x),1));
	mult_overflow61 = (zeros(length(x),1));
	mult_overflow62 = (zeros(length(x),1));
	mult_overflow63 = (zeros(length(x),1));
	mult_overflow64 = (zeros(length(x),1));
	mult_overflow65 = (zeros(length(x),1));
	mult_overflow66 = (zeros(length(x),1));
	mult_overflow67 = (zeros(length(x),1));
	mult_overflow68 = (zeros(length(x),1));
	mult_overflow69 = (zeros(length(x),1));
	mult_overflow70 = (zeros(length(x),1));
	mult_overflow71 = (zeros(length(x),1));
	mult_overflow72 = (zeros(length(x),1));

    sum1 = int32(zeros(length(x),1));
    sum2 = int32(zeros(length(x),1));
	sum3 = int32(zeros(length(x),1));
	sum4 = int32(zeros(length(x),1));
	sum5 = int32(zeros(length(x),1));
	sum6 = int32(zeros(length(x),1));
	sum7 = int32(zeros(length(x),1));
	sum8 = int32(zeros(length(x),1));
	sum9 = int32(zeros(length(x),1));
	sum10 = int32(zeros(length(x),1));
	sum11 = int32(zeros(length(x),1));
	sum12 = int32(zeros(length(x),1));
	sum13 = int32(zeros(length(x),1));
	sum14 = int32(zeros(length(x),1));
	sum15 = int32(zeros(length(x),1));
	sum16 = int32(zeros(length(x),1));
	sum17 = int32(zeros(length(x),1));
	sum18 = int32(zeros(length(x),1));
	sum19 = int32(zeros(length(x),1));
	sum20 = int32(zeros(length(x),1));
	sum21 = int32(zeros(length(x),1));
	sum22 = int32(zeros(length(x),1));
	sum23 = int32(zeros(length(x),1));
	sum24 = int32(zeros(length(x),1));
	sum25 = int32(zeros(length(x),1));
	sum26 = int32(zeros(length(x),1));
	sum27 = int32(zeros(length(x),1));
	sum28 = int32(zeros(length(x),1));
	sum29 = int32(zeros(length(x),1));
	sum30 = int32(zeros(length(x),1));
	sum31 = int32(zeros(length(x),1));
	sum32 = int32(zeros(length(x),1));
	sum33 = int32(zeros(length(x),1));
	sum34 = int32(zeros(length(x),1));
	sum35 = int32(zeros(length(x),1));
	sum36 = int32(zeros(length(x),1));
	sum37 = int32(zeros(length(x),1));
	sum38 = int32(zeros(length(x),1));
	sum39 = int32(zeros(length(x),1));
	sum40 = int32(zeros(length(x),1));
	sum41 = int32(zeros(length(x),1));
	sum42 = int32(zeros(length(x),1));
	sum43 = int32(zeros(length(x),1));
	sum44 = int32(zeros(length(x),1));
	sum45 = int32(zeros(length(x),1));
	sum46 = int32(zeros(length(x),1));
	sum47 = int32(zeros(length(x),1));
	sum48 = int32(zeros(length(x),1));
	sum49 = int32(zeros(length(x),1));
	sum50 = int32(zeros(length(x),1));
	sum51 = int32(zeros(length(x),1));
	sum52 = int32(zeros(length(x),1));
	sum53 = int32(zeros(length(x),1));
	sum54 = int32(zeros(length(x),1));
	sum55 = int32(zeros(length(x),1));
	sum56 = int32(zeros(length(x),1));
	sum57 = int32(zeros(length(x),1));
	sum58 = int32(zeros(length(x),1));
	sum59 = int32(zeros(length(x),1));
	sum60 = int32(zeros(length(x),1));
	sum61 = int32(zeros(length(x),1));
	sum62 = int32(zeros(length(x),1));
	sum63 = int32(zeros(length(x),1));
	sum64 = int32(zeros(length(x),1));
	sum65 = int32(zeros(length(x),1));
	sum66 = int32(zeros(length(x),1));
	sum67 = int32(zeros(length(x),1));
	sum68 = int32(zeros(length(x),1));
	sum69 = int32(zeros(length(x),1));
	sum70 = int32(zeros(length(x),1));
	
	y1_shift_round = int32(zeros(length(x),1));
	
	sum_overflow1 = (zeros(length(x),1));
    sum_overflow2 = (zeros(length(x),1));
	sum_overflow3 = (zeros(length(x),1));
	sum_overflow4 = (zeros(length(x),1));
	sum_overflow5 = (zeros(length(x),1));
	sum_overflow6 = (zeros(length(x),1));
	sum_overflow7 = (zeros(length(x),1));
	sum_overflow8 = (zeros(length(x),1));
	sum_overflow9 = (zeros(length(x),1));
	sum_overflow10 = (zeros(length(x),1));
	sum_overflow11 = (zeros(length(x),1));
	sum_overflow12 = (zeros(length(x),1));
	sum_overflow13 = (zeros(length(x),1));
	sum_overflow14 = (zeros(length(x),1));
	sum_overflow15 = (zeros(length(x),1));
	sum_overflow16 = (zeros(length(x),1));
	sum_overflow17 = (zeros(length(x),1));
	sum_overflow18 = (zeros(length(x),1));
	sum_overflow19 = (zeros(length(x),1));
	sum_overflow20 = (zeros(length(x),1));
	sum_overflow21 = (zeros(length(x),1));
	sum_overflow22 = (zeros(length(x),1));
	sum_overflow23 = (zeros(length(x),1));
	sum_overflow24 = (zeros(length(x),1));
	sum_overflow25 = (zeros(length(x),1));
	sum_overflow26 = (zeros(length(x),1));
	sum_overflow27 = (zeros(length(x),1));
	sum_overflow28 = (zeros(length(x),1));
	sum_overflow29 = (zeros(length(x),1));
	sum_overflow30 = (zeros(length(x),1));
	sum_overflow31 = (zeros(length(x),1));
	sum_overflow32 = (zeros(length(x),1));
	sum_overflow33 = (zeros(length(x),1));
	sum_overflow34 = (zeros(length(x),1));
	sum_overflow35 = (zeros(length(x),1));
	sum_overflow36 = (zeros(length(x),1));
	sum_overflow37 = (zeros(length(x),1));
	sum_overflow38 = (zeros(length(x),1));
	sum_overflow39 = (zeros(length(x),1));
	sum_overflow40 = (zeros(length(x),1));
	sum_overflow41 = (zeros(length(x),1));
	sum_overflow42 = (zeros(length(x),1));
	sum_overflow43 = (zeros(length(x),1));
	sum_overflow44 = (zeros(length(x),1));
	sum_overflow45 = (zeros(length(x),1));
	sum_overflow46 = (zeros(length(x),1));
	sum_overflow47 = (zeros(length(x),1));
	sum_overflow48 = (zeros(length(x),1));
	sum_overflow49 = (zeros(length(x),1));
	sum_overflow50 = (zeros(length(x),1));
	sum_overflow51 = (zeros(length(x),1));
	sum_overflow52 = (zeros(length(x),1));
	sum_overflow53 = (zeros(length(x),1));
	sum_overflow54 = (zeros(length(x),1));
	sum_overflow55 = (zeros(length(x),1));
	sum_overflow56 = (zeros(length(x),1));
	sum_overflow57 = (zeros(length(x),1));
	sum_overflow58 = (zeros(length(x),1));
	sum_overflow59 = (zeros(length(x),1));
	sum_overflow60 = (zeros(length(x),1));
	sum_overflow61 = (zeros(length(x),1));
	sum_overflow62 = (zeros(length(x),1));
	sum_overflow63 = (zeros(length(x),1));
	sum_overflow64 = (zeros(length(x),1));
	sum_overflow65 = (zeros(length(x),1));
	sum_overflow66 = (zeros(length(x),1));
	sum_overflow67 = (zeros(length(x),1));
	sum_overflow68 = (zeros(length(x),1));
	sum_overflow69 = (zeros(length(x),1));
	sum_overflow70 = (zeros(length(x),1));


    for n=1:nx
        p=p+1; if p>nb, p=1; end
        z(p) = x(n);
        acc = fi(0,1,width_out,0);
        k = p;
        for j=1:nb
            acc(:) = acc + b(j,:) * (z(k));
            k=k-1; if k<1, k=nb; end
        end

       y(n) = acc;

       %% integer
       %%
        buffer = [x(n) buffer(1:end-1)];

        [mult2(n), mult_overflow2(n)] = mult(coeff_int(2), int32(buffer(2)), 32);
        [mult3(n), mult_overflow3(n)] = mult(coeff_int(3), int32(buffer(3)), 32);
        [mult4(n), mult_overflow4(n)] = mult(coeff_int(4), int32(buffer(4)), 32);
        [mult5(n), mult_overflow5(n)] = mult(coeff_int(5), int32(buffer(5)), 32);
        [mult6(n), mult_overflow6(n)] = mult(coeff_int(6), int32(buffer(6)), 32);
        [mult7(n), mult_overflow7(n)] = mult(coeff_int(7), int32(buffer(7)), 32);
        [mult8(n), mult_overflow8(n)] = mult(coeff_int(8), int32(buffer(8)), 32);
        [mult9(n), mult_overflow9(n)] = mult(coeff_int(9), int32(buffer(9)), 32);
		
        [mult10(n), mult_overflow10(n)] = mult(coeff_int(10), int32(buffer(10)), 32);
        [mult11(n), mult_overflow11(n)] = mult(coeff_int(11), int32(buffer(11)), 32);
        [mult12(n), mult_overflow12(n)] = mult(coeff_int(12), int32(buffer(12)), 32);
        [mult13(n), mult_overflow13(n)] = mult(coeff_int(13), int32(buffer(13)), 32);
        [mult14(n), mult_overflow14(n)] = mult(coeff_int(14), int32(buffer(14)), 32);
        [mult15(n), mult_overflow15(n)] = mult(coeff_int(15), int32(buffer(15)), 32);
        [mult16(n), mult_overflow16(n)] = mult(coeff_int(16), int32(buffer(16)), 32);
        [mult17(n), mult_overflow17(n)] = mult(coeff_int(17), int32(buffer(17)), 32);
        [mult18(n), mult_overflow18(n)] = mult(coeff_int(18), int32(buffer(18)), 32);
        [mult19(n), mult_overflow19(n)] = mult(coeff_int(19), int32(buffer(19)), 32);

        [mult20(n), mult_overflow20(n)] = mult(coeff_int(20), int32(buffer(20)), 32);
        [mult21(n), mult_overflow21(n)] = mult(coeff_int(21), int32(buffer(21)), 32);
        [mult22(n), mult_overflow22(n)] = mult(coeff_int(22), int32(buffer(22)), 32);
        [mult23(n), mult_overflow23(n)] = mult(coeff_int(23), int32(buffer(23)), 32);
        [mult24(n), mult_overflow24(n)] = mult(coeff_int(24), int32(buffer(24)), 32);
        [mult25(n), mult_overflow25(n)] = mult(coeff_int(25), int32(buffer(25)), 32);
        [mult26(n), mult_overflow26(n)] = mult(coeff_int(26), int32(buffer(26)), 32);
        [mult27(n), mult_overflow27(n)] = mult(coeff_int(27), int32(buffer(27)), 32);
        [mult28(n), mult_overflow28(n)] = mult(coeff_int(28), int32(buffer(28)), 32);
        [mult29(n), mult_overflow29(n)] = mult(coeff_int(29), int32(buffer(29)), 32);

        [mult30(n), mult_overflow30(n)] = mult(coeff_int(30), int32(buffer(30)), 32);
        [mult31(n), mult_overflow31(n)] = mult(coeff_int(31), int32(buffer(31)), 32);
        [mult32(n), mult_overflow32(n)] = mult(coeff_int(32), int32(buffer(32)), 32);
        [mult33(n), mult_overflow33(n)] = mult(coeff_int(33), int32(buffer(33)), 32);
        [mult34(n), mult_overflow34(n)] = mult(coeff_int(34), int32(buffer(34)), 32);
        [mult35(n), mult_overflow35(n)] = mult(coeff_int(35), int32(buffer(35)), 32);
        [mult36(n), mult_overflow36(n)] = mult(coeff_int(36), int32(buffer(36)), 32);
        [mult37(n), mult_overflow37(n)] = mult(coeff_int(37), int32(buffer(37)), 32);
        [mult38(n), mult_overflow38(n)] = mult(coeff_int(38), int32(buffer(38)), 32);
        [mult39(n), mult_overflow39(n)] = mult(coeff_int(39), int32(buffer(39)), 32);

        [mult40(n), mult_overflow40(n)] = mult(coeff_int(40), int32(buffer(40)), 32);
        [mult41(n), mult_overflow41(n)] = mult(coeff_int(41), int32(buffer(41)), 32);
        [mult42(n), mult_overflow42(n)] = mult(coeff_int(42), int32(buffer(42)), 32);
        [mult43(n), mult_overflow43(n)] = mult(coeff_int(43), int32(buffer(43)), 32);
        [mult44(n), mult_overflow44(n)] = mult(coeff_int(44), int32(buffer(44)), 32);
        [mult45(n), mult_overflow45(n)] = mult(coeff_int(45), int32(buffer(45)), 32);
        [mult46(n), mult_overflow46(n)] = mult(coeff_int(46), int32(buffer(46)), 32);
        [mult47(n), mult_overflow47(n)] = mult(coeff_int(47), int32(buffer(47)), 32);
        [mult48(n), mult_overflow48(n)] = mult(coeff_int(48), int32(buffer(48)), 32);
        [mult49(n), mult_overflow49(n)] = mult(coeff_int(49), int32(buffer(49)), 32);

        [mult50(n), mult_overflow50(n)] = mult(coeff_int(50), int32(buffer(50)), 32);
        [mult51(n), mult_overflow51(n)] = mult(coeff_int(51), int32(buffer(51)), 32);
        [mult52(n), mult_overflow52(n)] = mult(coeff_int(52), int32(buffer(52)), 32);
        [mult53(n), mult_overflow53(n)] = mult(coeff_int(53), int32(buffer(53)), 32);
        [mult54(n), mult_overflow54(n)] = mult(coeff_int(54), int32(buffer(54)), 32);
        [mult55(n), mult_overflow55(n)] = mult(coeff_int(55), int32(buffer(55)), 32);
        [mult56(n), mult_overflow56(n)] = mult(coeff_int(56), int32(buffer(56)), 32);
        [mult57(n), mult_overflow57(n)] = mult(coeff_int(57), int32(buffer(57)), 32);
        [mult58(n), mult_overflow58(n)] = mult(coeff_int(58), int32(buffer(58)), 32);
        [mult59(n), mult_overflow59(n)] = mult(coeff_int(59), int32(buffer(59)), 32);

        [mult60(n), mult_overflow60(n)] = mult(coeff_int(60), int32(buffer(60)), 32);
        [mult61(n), mult_overflow61(n)] = mult(coeff_int(61), int32(buffer(61)), 32);
        [mult62(n), mult_overflow62(n)] = mult(coeff_int(62), int32(buffer(62)), 32);
        [mult63(n), mult_overflow63(n)] = mult(coeff_int(63), int32(buffer(63)), 32);
        [mult64(n), mult_overflow64(n)] = mult(coeff_int(64), int32(buffer(64)), 32);
        [mult65(n), mult_overflow65(n)] = mult(coeff_int(65), int32(buffer(65)), 32);
        [mult66(n), mult_overflow66(n)] = mult(coeff_int(66), int32(buffer(66)), 32);
        [mult67(n), mult_overflow67(n)] = mult(coeff_int(67), int32(buffer(67)), 32);
        [mult68(n), mult_overflow68(n)] = mult(coeff_int(68), int32(buffer(68)), 32);
        [mult69(n), mult_overflow69(n)] = mult(coeff_int(69), int32(buffer(69)), 32);

        [mult70(n), mult_overflow70(n)] = mult(coeff_int(70), int32(buffer(70)), 32);
        [mult71(n), mult_overflow71(n)] = mult(coeff_int(71), int32(buffer(71)), 32);
        [mult72(n), mult_overflow72(n)] = mult(coeff_int(72), int32(buffer(72)), 32);


        % mult2(n) = fi(2,0,2,0) * buffer(2); % fi(1,14,0)
        % mult2(n) = bitshift(mult2(n),shift);

        % mult3(n) = fi(-6,1,4,0) * buffer(3); % fi(1,16,0)
        % mult3(n) = bitshift(mult3(n),shift);

		% mult4(n) = fi(13,0,4,0) * buffer(4); % fi(1,16,0)
        % mult4(n) = bitshift(mult4(n),shift);

        % mult5(n) = fi(-24,1,6,0) * buffer(5); % fi(1,18,0)
        % mult5(n) = bitshift(mult5(n),shift);

        % mult6(n) = fi(39,0,6,0) * buffer(6); % fi(1,18,0)
        % mult6(n) = bitshift(mult6(n),shift);
		
		% mult7(n) = fi(-58,1,7,0) * buffer(7); % fi(1,19,0)
        % mult7(n) = bitshift(mult7(n),shift);
		
        % mult8(n) = fi(84,0,7,0) * buffer(8); % fi(1,19,0)
        % mult8(n) = bitshift(mult8(n),shift);
		
		% mult9(n) = fi(-116,1,8,0) * buffer(9); % fi(1,20,0)
        % mult9(n) = bitshift(mult9(n),shift);
		
		% mult10(n) = fi(156,0,8,0) * buffer(10); % fi(1,20,0)		
        % mult10(n) = bitshift(mult10(n),shift);
        
		% mult11(n) = fi(-206,1,9,0) * buffer(11); %  fi(1,21,0)		
        % mult11(n) = bitshift(mult11(n),shift);
		
		% mult12(n) = fi(266,0,9,0) * buffer(12); % fi (1,21,0)		
        % mult12(n) = bitshift(mult12(n),shift);
		
		% mult13(n) = fi(-339,1,10,0) * buffer(13); % fi(1,22,0)		
        % mult13(n) = bitshift(mult13(n),shift);
		
		% mult14(n) = fi(427,0,9,0) * buffer(14); % fi(1,22,0)		
        % mult14(n) = bitshift(mult14(n),shift);
		
		% mult15(n) = fi(-531,1,11,0) * buffer(15); % fi(1,23,0)		
        % mult15(n) = bitshift(mult15(n),shift);
		
		% mult16(n) = fi(655,0,10,0) * buffer(16); % fi(1,22,0)
        % mult16(n) = bitshift(mult16(n),shift);
		
		% mult17(n) = fi(-800,1,11,0) * buffer(17); % fi(1,23,0)
        % mult17(n) = bitshift(mult17(n),shift);
		
		% mult18(n) = fi(969,0,10,0) * buffer(18); % fi(1,22,0)
        % mult18(n) = bitshift(mult18(n),shift);
		
		% mult19(n) = fi(-1167,1,12,0) * buffer(19); % fi(1,24,0)
        % mult19(n) = bitshift(mult19(n),shift);
		
		% mult20(n) = fi(1396,0,11,0) * buffer(20); % fi(1,23,0)
		% mult20(n) = bitshift(mult20(n),shift);
		
		% mult21(n) = fi(-1662,1,12,0) * buffer(21); % fi(1,24,0)
		% mult21(n) = bitshift(mult21(n),shift);
		
		% mult22(n) = fi(1970,0,11,0) * buffer(22); % fi(1,23,0)
        % mult22(n) = bitshift(mult22(n),shift);
		
		% mult23(n) = fi(-2327,1,13,0) * buffer(23); % fi(1,25,0)
        % mult23(n) = bitshift(mult23(n),shift);
		
		% mult24(n) = fi(2742,0,12,0) * buffer(24); % fi(1,24,0)		
        % mult24(n) = bitshift(mult24(n),shift);
		
		% mult25(n) = fi(-3226,1,13,0) * buffer(25); % fi(1,25,0)		
        % mult25(n) = bitshift(mult25(n),shift);
		
		% mult26(n) = fi(3796,0,12,0) * buffer(26); % fi(1,24,0)		
        % mult26(n) = bitshift(mult26(n),shift);
		
		% mult27(n) = fi(-4474,1,14,0) * buffer(27); % fi(1,26,0)        
		% mult27(n) = bitshift(mult27(n),shift);
		
		% mult28(n) = fi(5291,0,13,0) * buffer(28); % fi(1,25,0)		
        % mult28(n) = bitshift(mult28(n),shift);
		
		% mult29(n) = fi(-6298,1,14,0) * buffer(29); % fi(1,26,0)		
        % mult29(n) = bitshift(mult29(n),shift);
		
		% mult30(n) = fi(7573,0,13,0) * buffer(30); % fi(1,25,0)		
        % mult30(n) = bitshift(mult30(n),shift);
		
		% mult31(n) = fi(-9249,1,15,0) * buffer(31); % fi(1,27,0)		
        % mult31(n) = bitshift(mult31(n),shift);
		
		% mult32(n) = fi(11573,0,14,0) * buffer(32); % fi(1,26,0)		
        % mult32(n) = bitshift(mult32(n),shift);
		
		% mult33(n) = fi(-15057,1,15,0) * buffer(33); % fi(1,27,0)		
        % mult33(n) = bitshift(mult33(n),shift);
		
		% mult34(n) = fi(20954,0,15,0) * buffer(34); % fi(1,27,0)		
        % mult34(n) = bitshift(mult34(n),shift);
		
		% mult35(n) = fi(-33395,1,17,0) * buffer(35); % fi(1,29,0)		
        % mult35(n) = bitshift(mult35(n),shift);
		
		% mult36(n) = fi(78533,0,17,0) * buffer(36); % fi(1,29,0)		
        % mult36(n) = bitshift(mult36(n),shift);
		
		% mult37(n) = fi(235966,0,18,0) * buffer(37); % fi(1,30,0)      
        %  mult37(n) = bitshift(mult37(n),shift);
		
		% mult38(n) = fi(-46973,1,17,0) * buffer(38); % fi(1,29,0)        
        % mult38(n) = bitshift(mult38(n),shift);
		
		% mult39(n) = fi(25812,0,15,0) * buffer(39); % fi(1,27,0)		
        % mult39(n) = bitshift(mult39(n),shift);
		
		% mult40(n) = fi(-17565,1,16,0) * buffer(40); % fi(1,28,0)        
        % mult40(n) = bitshift(mult40(n),shift);
		
		% mult41(n) = fi(13119,0,14,0) * buffer(41); % fi(1,26,0)        
        % mult41(n) = bitshift(mult41(n),shift);
		
		% mult42(n) = fi(-10307,1,15,0) * buffer(42); 		
        % mult42(n) = bitshift(mult42(n),shift);
		
		% mult43(n) = fi(8349,0,14,0) * buffer(43);        
        % mult43(n) = bitshift(mult43(n),shift);
		
		% mult44(n) = fi(-6895,1,14,0) * buffer(44);        
        % mult44(n) = bitshift(mult44(n),shift);
		
		% mult45(n) = fi(5767,0,13,0) * buffer(45);	        
        % mult45(n) = bitshift(mult45(n),shift);
		
		% mult46(n) = fi(-4862,1,14,0) * buffer(46);        
        % mult46(n) = bitshift(mult46(n),shift);
		
		% mult47(n) = fi(4120,0,13,0) * buffer(47);		
        % mult47(n) = bitshift(mult47(n),shift);
		
		% mult48(n) = fi(-3499,1,13,0) * buffer(48);       
        % mult48(n) = bitshift(mult48(n),shift);
		
		% mult49(n) = fi(2974,0,12,0) * buffer(49);       
        % mult49(n) = bitshift(mult49(n),shift);
		
		% mult50(n) = fi(-2526,1,13,0) * buffer(50);		
        % mult50(n) = bitshift(mult50(n),shift);
		
		% mult51(n) = fi(2142,0,12,0) * buffer(51);      
        % mult51(n) = bitshift(mult51(n),shift);
		
		% mult52(n) = fi(-1810,1,12,0) * buffer(52);      
        % mult52(n) = bitshift(mult52(n),shift);
		
		% mult53(n) = fi(1524,0,11,0) * buffer(53);		
        % mult53(n) = bitshift(mult53(n),shift);
		
		% mult54(n) = fi(-1277,1,12,0) * buffer(54);       
        % mult54(n) = bitshift(mult54(n),shift);
		
		% mult55(n) = fi(1064,0,11,0) * buffer(55);        
        % mult55(n) = bitshift(mult55(n),shift);
		
		% mult56(n) = fi(-881,1,11,0) * buffer(56);		
        % mult56(n) = bitshift(mult56(n),shift);
		
		% mult57(n) = fi(724,0,10,0) * buffer(57);        
        % mult57(n) = bitshift(mult57(n),shift);
		
		% mult58(n) = fi(-590,1,11,0) * buffer(58);        
        % mult58(n) = bitshift(mult58(n),shift);
		
		% mult59(n) = fi(477,0,9,0) * buffer(59);		
        % mult59(n) = bitshift(mult59(n),shift);
		
		% mult60(n) = fi(-381,1,10,0) * buffer(60);        
        % mult60(n) = bitshift(mult60(n),shift);
		
		% mult61(n) = fi(301,0,9,0) * buffer(61);        
        % mult61(n) = bitshift(mult61(n),shift);
		
		% mult62(n) = fi(-234,1,9,0) * buffer(62);		
        % mult62(n) = bitshift(mult62(n),shift);
		
		% mult63(n) = fi(180,0,8,0) * buffer(63);        
        % mult63(n) = bitshift(mult63(n),shift);
		
		% mult64(n) = fi(-135,1,9,0) * buffer(64);        
        % mult64(n) = bitshift(mult64(n),shift);
		
		% mult65(n) = fi(99,0,7,0) * buffer(65);		
        % mult65(n) = bitshift(mult65(n),shift);
		
		% mult66(n) = fi(-70,1,8,0) * buffer(66);		
        % mult66(n) = bitshift(mult66(n),shift);
		
		% mult67(n) = fi(48,0,6,0) * buffer(67);       
        % mult67(n) = bitshift(mult67(n),shift);
		
		% mult68(n) = fi(-31,1,6,0) * buffer(68);        
        % mult68(n) = bitshift(mult68(n),shift);
		
		% mult69(n) = fi(18,0,5,0) * buffer(69);		
        % mult69(n) = bitshift(mult69(n),shift);
		
		% mult70(n) = fi(-9,1,5,0) * buffer(70);        
        % mult70(n) = bitshift(mult70(n),shift);
		
		% mult71(n) = fi(4,0,3,0) * buffer(71);        
        % mult71(n) = bitshift(mult71(n),shift);
		
		% mult72(n) = fi(-1,1,2,0) * buffer(72);        
		% mult72(n) = bitshift(mult72(n),shift);
		

        %% adders
        %%

		% sum1(n) = mult2(n) + mult3(n);
		% sum1(n) = bitshift(sum1(n),shift_sum);

        [sum1(n), sum_overflow1(n)] = adder(mult2(n),  mult3(n), 32);
		
		[sum2(n), sum_overflow2(n)] = adder(sum1(n),  mult4(n), 32);
		[sum3(n), sum_overflow3(n)] = adder(sum2(n),  mult5(n), 32);
		[sum4(n), sum_overflow4(n)] = adder(sum3(n),  mult6(n), 32);
		[sum5(n), sum_overflow5(n)] = adder(sum4(n),  mult7(n), 32);
		[sum6(n), sum_overflow6(n)] = adder(sum5(n),  mult8(n), 32);
		[sum7(n), sum_overflow7(n)] = adder(sum6(n),  mult9(n), 32);
		[sum8(n), sum_overflow8(n)] = adder(sum7(n),  mult10(n), 32);
		[sum9(n), sum_overflow9(n)] = adder(sum8(n),  mult11(n), 32);
		
		[sum10(n), sum_overflow10(n)] = adder(sum9(n),  mult12(n), 32);
		[sum11(n), sum_overflow11(n)] = adder(sum10(n),  mult13(n), 32);
		[sum12(n), sum_overflow12(n)] = adder(sum11(n),  mult14(n), 32);
		[sum13(n), sum_overflow13(n)] = adder(sum12(n),  mult15(n), 32);
		[sum14(n), sum_overflow14(n)] = adder(sum13(n),  mult16(n), 32);
		[sum15(n), sum_overflow15(n)] = adder(sum14(n),  mult17(n), 32);
		[sum16(n), sum_overflow16(n)] = adder(sum15(n),  mult18(n), 32);
		[sum17(n), sum_overflow17(n)] = adder(sum16(n),  mult19(n), 32);
		[sum18(n), sum_overflow18(n)] = adder(sum17(n),  mult20(n), 32);
		[sum19(n), sum_overflow19(n)] = adder(sum18(n),  mult21(n), 32);
		
		[sum20(n), sum_overflow20(n)] = adder(sum19(n),  mult22(n), 32);
		[sum21(n), sum_overflow21(n)] = adder(sum20(n),  mult23(n), 32);
		[sum22(n), sum_overflow22(n)] = adder(sum21(n),  mult24(n), 32);
		[sum23(n), sum_overflow23(n)] = adder(sum22(n),  mult25(n), 32);
		[sum24(n), sum_overflow24(n)] = adder(sum23(n),  mult26(n), 32);
		[sum25(n), sum_overflow25(n)] = adder(sum24(n),  mult27(n), 32);
		[sum26(n), sum_overflow26(n)] = adder(sum25(n),  mult28(n), 32);
		[sum27(n), sum_overflow27(n)] = adder(sum26(n),  mult29(n), 32);
		[sum28(n), sum_overflow28(n)] = adder(sum27(n),  mult30(n), 32);
		[sum29(n), sum_overflow29(n)] = adder(sum28(n),  mult31(n), 32);
					
		[sum30(n), sum_overflow30(n)] = adder(sum29(n),  mult32(n), 32);
		[sum31(n), sum_overflow31(n)] = adder(sum30(n),  mult33(n), 32);
		[sum32(n), sum_overflow32(n)] = adder(sum31(n),  mult34(n), 32);
		[sum33(n), sum_overflow33(n)] = adder(sum32(n),  mult35(n), 32);
		[sum34(n), sum_overflow34(n)] = adder(sum33(n),  mult36(n), 32);
		[sum35(n), sum_overflow35(n)] = adder(sum34(n),  mult37(n), 32);
		[sum36(n), sum_overflow36(n)] = adder(sum35(n),  mult38(n), 32);
		[sum37(n), sum_overflow37(n)] = adder(sum36(n),  mult39(n), 32);
		[sum38(n), sum_overflow38(n)] = adder(sum37(n),  mult40(n), 32);
		[sum39(n), sum_overflow39(n)] = adder(sum38(n),  mult41(n), 32);
		
		[sum40(n), sum_overflow40(n)] = adder(sum39(n),  mult42(n), 32);
		[sum41(n), sum_overflow41(n)] = adder(sum40(n),  mult43(n), 32);
		[sum42(n), sum_overflow42(n)] = adder(sum41(n),  mult44(n), 32);
		[sum43(n), sum_overflow43(n)] = adder(sum42(n),  mult45(n), 32);
		[sum44(n), sum_overflow44(n)] = adder(sum43(n),  mult46(n), 32);
		[sum45(n), sum_overflow45(n)] = adder(sum44(n),  mult47(n), 32);
		[sum46(n), sum_overflow46(n)] = adder(sum45(n),  mult48(n), 32);
		[sum47(n), sum_overflow47(n)] = adder(sum46(n),  mult49(n), 32);
		[sum48(n), sum_overflow48(n)] = adder(sum47(n),  mult50(n), 32);
		[sum49(n), sum_overflow49(n)] = adder(sum48(n),  mult51(n), 32);	
		
		[sum50(n), sum_overflow50(n)] = adder(sum49(n),  mult52(n), 32);
		[sum51(n), sum_overflow51(n)] = adder(sum50(n),  mult53(n), 32);
		[sum52(n), sum_overflow52(n)] = adder(sum51(n),  mult54(n), 32);
		[sum53(n), sum_overflow53(n)] = adder(sum52(n),  mult55(n), 32);
		[sum54(n), sum_overflow54(n)] = adder(sum53(n),  mult56(n), 32);
		[sum55(n), sum_overflow55(n)] = adder(sum54(n),  mult57(n), 32);
		[sum56(n), sum_overflow56(n)] = adder(sum55(n),  mult58(n), 32);
		[sum57(n), sum_overflow57(n)] = adder(sum56(n),  mult59(n), 32);
		[sum58(n), sum_overflow58(n)] = adder(sum57(n),  mult60(n), 32);
		[sum59(n), sum_overflow59(n)] = adder(sum58(n),  mult61(n), 32);
		
		[sum60(n), sum_overflow60(n)] = adder(sum59(n),  mult62(n), 32);
		[sum61(n), sum_overflow61(n)] = adder(sum60(n),  mult63(n), 32);
		[sum62(n), sum_overflow62(n)] = adder(sum61(n),  mult64(n), 32);
		[sum63(n), sum_overflow63(n)] = adder(sum62(n),  mult65(n), 32);
		[sum64(n), sum_overflow64(n)] = adder(sum63(n),  mult66(n), 32);
		[sum65(n), sum_overflow65(n)] = adder(sum64(n),  mult67(n), 32);
		[sum66(n), sum_overflow66(n)] = adder(sum65(n),  mult68(n), 32);
		[sum67(n), sum_overflow67(n)] = adder(sum66(n),  mult69(n), 32);
		[sum68(n), sum_overflow68(n)] = adder(sum67(n),  mult70(n), 32);
		[sum69(n), sum_overflow69(n)] = adder(sum68(n),  mult71(n), 32);
		[sum70(n), sum_overflow70(n)] = adder(sum69(n),  mult72(n), 32);

        y1_shift_round(n) = bitshift(sum70(n),-13);
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		% sum2(n) = sum1(n) + mult4(n);	
		% sum2(n) = bitshift(sum2(n),shift_sum);
		
		% sum3(n) = sum2(n) + mult5(n);
		% sum3(n) = bitshift(sum3(n),shift_sum);

		% sum4(n) = sum3(n) + mult6(n);
		% sum4(n) = bitshift(sum4(n),shift_sum);

		% sum5(n) = sum4(n) + mult7(n);
		% sum5(n) = bitshift(sum5(n),shift_sum);
		
        % sum6(n) = sum5(n) + mult8(n);
		% sum6(n) = bitshift(sum6(n),shift_sum);
		
        % sum7(n) = sum6(n) + mult9(n);
		% sum7(n) = bitshift(sum7(n),shift_sum);
		
        % sum8(n) = sum7(n) + mult10(n);	
        % sum8(n) = bitshift(sum8(n),shift_sum);
		
        % sum9(n) = sum8(n) + mult11(n);
        % sum9(n) = bitshift(sum9(n),shift_sum);
		
        % sum10(n) = sum9(n) + mult12(n);
		% sum10(n) = bitshift(sum10(n),shift_sum);
		
        % sum11(n) = sum10(n) + mult13(n);
		% sum11(n) = bitshift(sum11(n),shift_sum);
		
        % sum12(n) = sum11(n) + mult14(n);
		% sum12(n) = bitshift(sum12(n),shift_sum);
		
        % sum13(n) = sum12(n) + mult15(n);
		% sum13(n) = bitshift(sum13(n),shift_sum);
		
        % sum14(n) = sum13(n) + mult16(n);
		% sum14(n) = bitshift(sum14(n),shift_sum);
		
        % sum15(n) = sum14(n) + mult17(n);
		% sum15(n) = bitshift(sum15(n),shift_sum);
		
        % sum16(n) = sum15(n) + mult18(n);	
        % sum16(n) = bitshift(sum16(n),shift_sum);
		
        % sum17(n) = sum16(n) + mult19(n);
        % sum17(n) = bitshift(sum17(n),shift_sum);
		
        % sum18(n) = sum17(n) + mult20(n);	
		% sum18(n) = bitshift(sum18(n),shift_sum);
		
        % sum19(n) = sum18(n) + mult21(n);
		% sum19(n) = bitshift(sum19(n),shift_sum);
		
        % sum20(n) = sum19(n) + mult22(n);	
        % sum20(n) = bitshift(sum20(n),shift_sum);
		
        % sum21(n) = sum20(n) + mult23(n);
		% sum21(n) = bitshift(sum21(n),shift_sum);
		
        % sum22(n) = sum21(n) + mult24(n);	
        % sum22(n) = bitshift(sum22(n),shift_sum);
		
        % sum23(n) = sum22(n) + mult25(n);
        % sum23(n) = bitshift(sum23(n),shift_sum);
		
        % sum24(n) = sum23(n) + mult26(n);	
        % sum24(n) = bitshift(sum24(n),shift_sum);

        % sum25(n) = sum24(n) + mult27(n);
		% sum25(n) = bitshift(sum25(n),shift_sum);
		
        % sum26(n) = sum25(n) + mult28(n);	
        % sum26(n) = bitshift(sum26(n),shift_sum);
		
        % sum27(n) = sum26(n) + mult29(n);
        % sum27(n) = bitshift(sum27(n),shift_sum);
		
        % sum28(n) = sum27(n) + mult30(n);	
        % sum28(n) = bitshift(sum28(n),shift_sum);
		
        % sum29(n) = sum28(n) + mult31(n);
		% sum29(n) = bitshift(sum29(n),shift_sum);
		
        % sum30(n) = sum29(n) + mult32(n);	
        % sum30(n) = bitshift(sum30(n),shift_sum);
		
        % sum31(n) = sum30(n) + mult33(n);
        % sum31(n) = bitshift(sum31(n),shift_sum);

        % sum32(n) = sum31(n) + mult34(n);
		% sum32(n) = bitshift(sum32(n),shift_sum);
		
        % sum33(n) = sum32(n) + mult35(n);
        % sum33(n) = bitshift(sum33(n),shift_sum);
		
        % sum34(n) = sum33(n) + mult36(n);
		% sum34(n) = bitshift(sum34(n),shift_sum);
		
        % sum35(n) = sum34(n) + mult37(n);
		% sum35(n) = bitshift(sum35(n),shift_sum);
		
        % sum36(n) = sum35(n) + mult38(n);	
		% sum36(n) = bitshift(sum36(n),shift_sum);
		
        % sum37(n) = sum36(n) + mult39(n);
		% sum37(n) = bitshift(sum37(n),shift_sum);
		
        % sum38(n) = sum37(n) + mult40(n);	
		% sum38(n) = bitshift(sum38(n),shift_sum);
		
        % sum39(n) = sum38(n) + mult41(n);
		% sum39(n) = bitshift(sum39(n),shift_sum);
		
        % sum40(n) = sum39(n) + mult42(n);
        % sum40(n) = bitshift(sum40(n),shift_sum);
		
        % sum41(n) = sum40(n) + mult43(n);
        % sum41(n) = bitshift(sum41(n),shift_sum);
		
        % sum42(n) = sum41(n) + mult44(n);	
		% sum42(n) = bitshift(sum42(n),shift_sum);
		
        % sum43(n) = sum42(n) + mult45(n);
        % sum43(n) = bitshift(sum43(n),shift_sum);
		
        % sum44(n) = sum43(n) + mult46(n);	
        % sum44(n) = bitshift(sum44(n),shift_sum);
		
        % sum45(n) = sum44(n) + mult47(n);
        % sum45(n) = bitshift(sum45(n),shift_sum);
		
        % sum46(n) = sum45(n) + mult48(n);	
		% sum46(n) = bitshift(sum46(n),shift_sum);
		
        % sum47(n) = sum46(n) + mult49(n);
        % sum47(n) = bitshift(sum47(n),shift_sum);
		
        % sum48(n) = sum47(n) + mult50(n);	
        % sum48(n) = bitshift(sum48(n),shift_sum);
		
        % sum49(n) = sum48(n) + mult51(n);
        % sum49(n) = bitshift(sum49(n),shift_sum);
		
        % sum50(n) = sum49(n) + mult52(n);	
		% sum50(n) = bitshift(sum50(n),shift_sum);
		
        % sum51(n) = sum50(n) + mult53(n);
        % sum51(n) = bitshift(sum51(n),shift_sum);
		
        % sum52(n) = sum51(n) + mult54(n);	
        % sum52(n) = bitshift(sum52(n),shift_sum);
		
        % sum53(n) = sum52(n) + mult55(n);
        % sum53(n) = bitshift(sum53(n),shift_sum);
		
        % sum54(n) = sum53(n) + mult56(n);
		% sum54(n) = bitshift(sum54(n),shift_sum);
		
        % sum55(n) = sum54(n) + mult57(n);
		% sum55(n) = bitshift(sum55(n),shift_sum);
		
        % sum56(n) = sum55(n) + mult58(n);	
        % sum56(n) = bitshift(sum56(n),shift_sum);

        % sum57(n) = sum56(n) + mult59(n);
		% sum57(n) = bitshift(sum57(n),shift_sum);
		
        % sum58(n) = sum57(n) + mult60(n);	
        % sum58(n) = bitshift(sum58(n),shift_sum);
		
        % sum59(n) = sum58(n) + mult61(n);
        % sum59(n) = bitshift(sum59(n),shift_sum);
		
        % sum60(n) = sum59(n) + mult62(n);	
        % sum60(n) = bitshift(sum60(n),shift_sum);
		
        % sum61(n) = sum60(n) + mult63(n);
		% sum61(n) = bitshift(sum61(n),shift_sum);
		
        % sum62(n) = sum61(n) + mult64(n);	
        % sum62(n) = bitshift(sum62(n),shift_sum);
		
        % sum63(n) = sum62(n) + mult65(n);
        % sum63(n) = bitshift(sum63(n),shift_sum);
		
        % sum64(n) = sum63(n) + mult66(n);	
        % sum64(n) = bitshift(sum64(n),shift_sum);
		
        % sum65(n) = sum64(n) + mult67(n);
        % sum65(n) = bitshift(sum65(n),shift_sum);
		
        % sum66(n) = sum65(n) + mult68(n);	
        % sum66(n) = bitshift(sum66(n),shift_sum);
		
        % sum67(n) = sum66(n) + mult69(n);
        % sum67(n) = bitshift(sum67(n),shift_sum);
		
        % sum68(n) = sum67(n) + mult70(n);	
        % sum68(n) = bitshift(sum68(n),shift_sum);
		
        % sum69(n) = sum68(n) + mult71(n);
        % sum69(n) = bitshift(sum69(n),shift_sum);
		
        % y1(n) = sum69(n) + mult72(n);
        % y1_shift_round(n) = bitshift(y1(n),-13);
        
    end

    max_mult2 = max(mult2);
	max_mult3 = max(mult3);
	max_mult4 = max(mult4);
	max_mult5 = max(mult5);
	max_mult6 = max(mult6);
	max_mult7 = max(mult7);
	max_mult8 = max(mult8);
	max_mult9 = max(mult9);
	max_mult10 = max(mult10);
	max_mult11 = max(mult11);
	max_mult12 = max(mult12);
	max_mult13 = max(mult13);
	max_mult14 = max(mult14);
	max_mult15 = max(mult15);
	max_mult16 = max(mult16);
	max_mult17 = max(mult17);
	max_mult18 = max(mult18);
	max_mult19 = max(mult19);
	max_mult20 = max(mult20);
	max_mult21 = max(mult21);
	max_mult22 = max(mult22);
	max_mult23 = max(mult23);
	max_mult24 = max(mult24);
	max_mult25 = max(mult25);
	max_mult26 = max(mult26);
	max_mult27 = max(mult27);
	max_mult28 = max(mult28);
	max_mult29 = max(mult29);
	max_mult30 = max(mult30);
	max_mult31 = max(mult31);
	max_mult32 = max(mult32);
	max_mult33 = max(mult33);
	max_mult34 = max(mult34);
	max_mult35 = max(mult35);
	max_mult36 = max(mult36);
	max_mult37 = max(mult37);
	max_mult38 = max(mult38);
	max_mult39 = max(mult39);
	max_mult40 = max(mult40);
	max_mult41 = max(mult41);
	max_mult42 = max(mult42);
	max_mult43 = max(mult43);
	max_mult44 = max(mult44);
	max_mult45 = max(mult45);
	max_mult46 = max(mult46);
	max_mult47 = max(mult47);
	max_mult48 = max(mult48);
	max_mult49 = max(mult49);
	max_mult50 = max(mult50);
	max_mult51 = max(mult51);
	max_mult52 = max(mult52);
	max_mult53 = max(mult53);
	max_mult54 = max(mult54);
	max_mult55 = max(mult55);
	max_mult56 = max(mult56);
	max_mult57 = max(mult57);
	max_mult58 = max(mult58);
	max_mult59 = max(mult59);
	max_mult60 = max(mult60);
	max_mult61 = max(mult61);
	max_mult62 = max(mult62);
	max_mult63 = max(mult63);
	max_mult64 = max(mult64);
	max_mult65 = max(mult65);
	max_mult66 = max(mult66);
	max_mult67 = max(mult67);
	max_mult68 = max(mult68);
	max_mult69 = max(mult69);
	max_mult70 = max(mult70);
	max_mult71 = max(mult71);
	max_mult72 = max(mult72);
	
	max_sum1 = max(sum1);
	max_sum2 = max(sum2);
	max_sum3 = max(sum3);
	max_sum4 = max(sum4);
	max_sum5 = max(sum5);
	max_sum6 = max(sum6);
	max_sum7 = max(sum7);
	max_sum8 = max(sum8);
	max_sum9 = max(sum9);
	max_sum10 = max(sum10);
	max_sum11 = max(sum11);
	max_sum12 = max(sum12);
	max_sum13 = max(sum13);
	max_sum14 = max(sum14);
	max_sum15 = max(sum15);
	max_sum16 = max(sum16);
	max_sum17 = max(sum17);
	max_sum18 = max(sum18);
	max_sum19 = max(sum19);
	max_sum20 = max(sum20);
	max_sum21 = max(sum21);
	max_sum22 = max(sum22);
	max_sum23 = max(sum23);
	max_sum24 = max(sum24);
	max_sum25 = max(sum25);
	max_sum26 = max(sum26);
	max_sum27 = max(sum27);
	max_sum28 = max(sum28);
	max_sum29 = max(sum29);
	max_sum30 = max(sum30);
	max_sum31 = max(sum31);
	max_sum32 = max(sum32);
	max_sum33 = max(sum33);
	max_sum34 = max(sum34);
	max_sum35 = max(sum35);
	max_sum36 = max(sum36);
	max_sum37 = max(sum37);
	max_sum38 = max(sum38);
	max_sum39 = max(sum39);
	max_sum40 = max(sum40);
	max_sum41 = max(sum41);
	max_sum42 = max(sum42);
	max_sum43 = max(sum43);
	max_sum44 = max(sum44);
	max_sum45 = max(sum45);
	max_sum46 = max(sum46);
	max_sum47 = max(sum47);
	max_sum48 = max(sum48);
	max_sum49 = max(sum49);
	max_sum50 = max(sum50);
	max_sum51 = max(sum51);
	max_sum52 = max(sum52);
	max_sum53 = max(sum53);
	max_sum54 = max(sum54);
	max_sum55 = max(sum55);
	max_sum56 = max(sum56);
	max_sum57 = max(sum57);
	max_sum58 = max(sum58);
	max_sum59 = max(sum59);
	max_sum60 = max(sum60);
	max_sum61 = max(sum61);
	max_sum62 = max(sum62);
	max_sum63 = max(sum63);
	max_sum64 = max(sum64);
	max_sum65 = max(sum65);
	max_sum66 = max(sum66);
	max_sum67 = max(sum67);
	max_sum68 = max(sum68);
	max_sum69 = max(sum69);
	max_sum70 = max(sum70);

	
end
