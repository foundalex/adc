function [y, y1] = fir_filter(b, x, width_in, width_out)

    z = fi(zeros(size(b)),1,width_in,0);
    y = fi(zeros(size(x)),1,width_out,0);

    p = 0;
    nx = length(x);
    nb = length(b);

    buffer = zeros(1,73);

    shift = -18;

    for n=1:nx
        p=p+1; if p>nb, p=1; end
        z(p) = x(n);
        acc = fi(0,1,width_out,0);
        k = p;
        for j=1:nb
            acc(:) = acc + b(j,:) * (z(k));
            k=k-1; if k<1, k=nb; end
        end

        buffer = [x(n) buffer(1:end-1)];

        mult2(n) = fi(2,0,2,0) * buffer(2); % max 4094
        mult2(n) = bitshift(mult2(n),shift);

        mult3(n) = fi(-6,1,4,0) * buffer(3); % max 12288
        mult3(n) = bitshift(mult3(n),shift);

		mult4(n) = fi(13,0,4,0) * buffer(4); % max 26611
        mult4(n) = bitshift(mult4(n),shift);

        mult5(n) = fi(-24,1,6,0) * buffer(5);
        mult5(n) = bitshift(mult5(n),shift);

        mult6(n) = fi(39,0,6,0) * buffer(6);
		mult6(n) = bitshift(mult6(n),shift);
		
		mult7(n) = fi(-58,1,7,0) * buffer(7);
		mult7(n) = bitshift(mult7(n),shift);
		
        mult8(n) = fi(84,0,7,0) * buffer(8);
        mult8(n) = bitshift(mult8(n),shift);
		
		mult9(n) = fi(-116,1,8,0) * buffer(9);
		mult9(n) = bitshift(mult9(n),shift);
		
		mult10(n) = fi(156,0,8,0) * buffer(10);
		mult10(n) = bitshift(mult10(n),shift);
        
		mult11(n) = fi(-206,1,9,0) * buffer(11);
		mult11(n) = bitshift(mult11(n),shift);
		
		mult12(n) = fi(266,0,9,0) * buffer(12);
		mult12(n) = bitshift(mult12(n),shift);
		
		mult13(n) = fi(-339,1,10,0) * buffer(13);
		mult13(n) = bitshift(mult13(n),shift);
		
		mult14(n) = fi(427,0,9,0) * buffer(14);
		mult14(n) = bitshift(mult14(n),shift);
		
		mult15(n) = fi(-531,1,11,0) * buffer(15);
		mult15(n) = bitshift(mult15(n),shift);
		
		mult16(n) = fi(655,0,10,0) * buffer(16);
		mult16(n) = bitshift(mult16(n),shift);
		
		mult17(n) = fi(-800,1,11,0) * buffer(17);
		mult17(n) = bitshift(mult17(n),shift);
		
		mult18(n) = fi(969,0,10,0) * buffer(18);
		mult18(n) = bitshift(mult18(n),shift);
		
		mult19(n) = fi(-1167,1,12,0) * buffer(19);
		mult19(n) = bitshift(mult19(n),shift);
		
		mult20(n) = fi(1396,0,11,0) * buffer(20);
		mult20(n) = bitshift(mult20(n),shift);
		
		mult21(n) = fi(-1662,1,12,0) * buffer(21);
		mult21(n) = bitshift(mult21(n),shift);
		
		mult22(n) = fi(1970,0,11,0) * buffer(22);
		mult22(n) = bitshift(mult22(n),shift);
		
		mult23(n) = fi(-2327,1,13,0) * buffer(23);
		mult23(n) = bitshift(mult23(n),shift);
		
		mult24(n) = fi(2742,0,12,0) * buffer(24);
		mult24(n) = bitshift(mult24(n),shift);
		
		mult25(n) = fi(-3226,1,13,0) * buffer(25);
		mult25(n) = bitshift(mult25(n),shift);
		
		mult26(n) = fi(3796,0,12,0) * buffer(26);
		mult26(n) = bitshift(mult26(n),shift);
		
		mult27(n) = fi(-4474,1,14,0) * buffer(27);
		mult27(n) = bitshift(mult27(n),shift);
		
		mult28(n) = fi(5291,0,13,0) * buffer(28);
		mult28(n) = bitshift(mult28(n),shift);
		
		mult29(n) = fi(-6298,1,14,0) * buffer(29);
		mult29(n) = bitshift(mult29(n),shift);
		
		mult30(n) = fi(7573,0,13,0) * buffer(30);
		mult30(n) = bitshift(mult30(n),shift);
		
		mult31(n) = fi(-9249,1,15,0) * buffer(31);
		mult31(n) = bitshift(mult31(n),shift);
		
		mult32(n) = fi(11573,0,14,0) * buffer(32);
		mult32(n) = bitshift(mult32(n),shift);
		
		mult33(n) = fi(-15057,1,15,0) * buffer(33);
		mult33(n) = bitshift(mult33(n),shift);
		
		mult34(n) = fi(20954,0,15,0) * buffer(34);
		mult34(n) = bitshift(mult34(n),shift);
		
		mult35(n) = fi(-33395,1,17,0) * buffer(35);
		mult35(n) = bitshift(mult35(n),shift);
		
		mult36(n) = fi(78533,0,17,0) * buffer(36);
		mult36(n) = bitshift(mult36(n),shift);
		
		mult37(n) = fi(235966,0,18,0) * buffer(37);
        mult37(n) = bitshift(mult37(n),shift);
		
		mult38(n) = fi(-46973,1,17,0) * buffer(38);
        mult38(n) = bitshift(mult38(n),shift);
		
		mult39(n) = fi(25812,0,15,0) * buffer(39);
		mult39(n) = bitshift(mult39(n),shift);
		
		mult40(n) = fi(-17565,1,16,0) * buffer(40);
        mult40(n) = bitshift(mult40(n),shift);
		
		mult41(n) = fi(13119,0,14,0) * buffer(41);
        mult41(n) = bitshift(mult41(n),shift);
		
		mult42(n) = fi(-10307,1,15,0) * buffer(42);
		mult42(n) = bitshift(mult42(n),shift);
		
		mult43(n) = fi(8349,0,14,0) * buffer(43);
        mult43(n) = bitshift(mult43(n),shift);
		
		mult44(n) = fi(-6895,1,14,0) * buffer(44);
        mult44(n) = bitshift(mult44(n),shift);
		
		mult45(n) = fi(5767,0,13,0) * buffer(45);	
        mult45(n) = bitshift(mult45(n),shift);
		
		mult46(n) = fi(-4862,1,14,0) * buffer(46);
        mult46(n) = bitshift(mult46(n),shift);
		
		mult47(n) = fi(4120,0,13,0) * buffer(47);
		mult47(n) = bitshift(mult47(n),shift);
		
		mult48(n) = fi(-3499,1,13,0) * buffer(48);
        mult48(n) = bitshift(mult48(n),shift);
		
		mult49(n) = fi(2974,0,12,0) * buffer(49);
        mult49(n) = bitshift(mult49(n),shift);
		
		mult50(n) = fi(-2526,1,13,0) * buffer(50);
		mult50(n) = bitshift(mult50(n),shift);
		
		mult51(n) = fi(2142,0,12,0) * buffer(51);
        mult51(n) = bitshift(mult51(n),shift);
		
		mult52(n) = fi(-1810,1,12,0) * buffer(52);
        mult52(n) = bitshift(mult52(n),shift);
		
		mult53(n) = fi(1524,0,11,0) * buffer(53);
		mult53(n) = bitshift(mult53(n),shift);
		
		mult54(n) = fi(-1277,1,12,0) * buffer(54);
        mult54(n) = bitshift(mult54(n),shift);
		
		mult55(n) = fi(1064,0,11,0) * buffer(55);
        mult55(n) = bitshift(mult55(n),shift);
		
		mult56(n) = fi(-881,1,11,0) * buffer(56);
		mult56(n) = bitshift(mult56(n),shift);
		
		mult57(n) = fi(724,0,10,0) * buffer(57);
        mult57(n) = bitshift(mult57(n),shift);
		
		mult58(n) = fi(-590,1,11,0) * buffer(58);
        mult58(n) = bitshift(mult58(n),shift);
		
		mult59(n) = fi(477,0,9,0) * buffer(59);
		mult59(n) = bitshift(mult59(n),shift);
		
		mult60(n) = fi(-381,1,10,0) * buffer(60);
        mult60(n) = bitshift(mult60(n),shift);
		
		mult61(n) = fi(301,0,9,0) * buffer(61);
        mult61(n) = bitshift(mult61(n),shift);
		
		mult62(n) = fi(-234,1,9,0) * buffer(62);
		mult62(n) = bitshift(mult62(n),shift);
		
		mult63(n) = fi(180,0,8,0) * buffer(63);
        mult63(n) = bitshift(mult63(n),shift);
		
		mult64(n) = fi(-135,1,9,0) * buffer(64);
        mult64(n) = bitshift(mult64(n),shift);
		
		mult65(n) = fi(99,0,7,0) * buffer(65);
		mult65(n) = bitshift(mult65(n),shift);
		
		mult66(n) = fi(-70,1,8,0) * buffer(66);
		mult66(n) = bitshift(mult66(n),shift);
		
		mult67(n) = fi(48,0,6,0) * buffer(67);
        mult67(n) = bitshift(mult67(n),shift);
		
		mult68(n) = fi(-31,1,6,0) * buffer(68);
        mult68(n) = bitshift(mult68(n),shift);
		
		mult69(n) = fi(18,0,5,0) * buffer(69);
		mult69(n) = bitshift(mult69(n),shift);
		
		mult70(n) = fi(-9,1,5,0) * buffer(70);
        mult70(n) = bitshift(mult70(n),shift);
		
		mult71(n) = fi(4,0,3,0) * buffer(71);
        mult71(n) = bitshift(mult71(n),shift);
		
		mult72(n) = fi(-1,1,2,0) * buffer(72);
		mult72(n) = bitshift(mult72(n),shift);
		
		
		
		
		
		sum1(n) = mult2(n) + mult3(n);
		sum2(n) = sum1(n) + mult4(n);	
		sum3(n) = sum2(n) + mult5(n);
		sum4(n) = sum3(n) + mult6(n);	
		sum5(n) = sum4(n) + mult7(n);
		sum6(n) = sum5(n) + mult8(n);	
		sum7(n) = sum6(n) + mult9(n);
		sum8(n) = sum7(n) + mult10(n);	
		sum9(n) = sum8(n) + mult11(n);
		sum10(n) = sum9(n) + mult12(n);	
		sum11(n) = sum10(n) + mult13(n);
		sum12(n) = sum11(n) + mult14(n);	
		sum13(n) = sum12(n) + mult15(n);
		sum14(n) = sum13(n) + mult16(n);	
		sum15(n) = sum14(n) + mult17(n);
		sum16(n) = sum15(n) + mult18(n);	
		sum17(n) = sum16(n) + mult19(n);
		sum18(n) = sum17(n) + mult20(n);	
		sum19(n) = sum18(n) + mult21(n);
		sum20(n) = sum19(n) + mult22(n);	
		sum21(n) = sum20(n) + mult23(n);
		sum22(n) = sum21(n) + mult24(n);	
		sum23(n) = sum22(n) + mult25(n);
		sum24(n) = sum23(n) + mult26(n);	
		sum25(n) = sum24(n) + mult27(n);
		sum26(n) = sum25(n) + mult28(n);	
		sum27(n) = sum26(n) + mult29(n);
		sum28(n) = sum27(n) + mult30(n);	
		sum29(n) = sum28(n) + mult31(n);
		sum30(n) = sum29(n) + mult32(n);	
		sum31(n) = sum30(n) + mult33(n);
		sum32(n) = sum31(n) + mult34(n);	
		sum33(n) = sum32(n) + mult35(n);
		sum34(n) = sum33(n) + mult36(n);	
		sum35(n) = sum34(n) + mult37(n);
		sum36(n) = sum35(n) + mult38(n);	
		sum37(n) = sum36(n) + mult39(n);
		sum38(n) = sum37(n) + mult40(n);	
		sum39(n) = sum38(n) + mult41(n);
		sum40(n) = sum39(n) + mult42(n);	
		sum41(n) = sum40(n) + mult43(n);
		sum42(n) = sum41(n) + mult44(n);	
		sum43(n) = sum42(n) + mult45(n);
		sum44(n) = sum43(n) + mult46(n);	
		sum45(n) = sum44(n) + mult47(n);
		sum46(n) = sum45(n) + mult48(n);	
		sum47(n) = sum46(n) + mult49(n);
		sum48(n) = sum47(n) + mult50(n);	
		sum49(n) = sum48(n) + mult51(n);
		sum50(n) = sum49(n) + mult52(n);	
		sum51(n) = sum50(n) + mult53(n);
		sum52(n) = sum51(n) + mult54(n);	
		sum53(n) = sum52(n) + mult55(n);
		sum54(n) = sum53(n) + mult56(n);	
		sum55(n) = sum54(n) + mult57(n);
		sum56(n) = sum55(n) + mult58(n);	
		sum57(n) = sum56(n) + mult59(n);
		sum58(n) = sum57(n) + mult60(n);	
		sum59(n) = sum58(n) + mult61(n);
		sum60(n) = sum59(n) + mult62(n);	
		sum61(n) = sum60(n) + mult63(n);
		sum62(n) = sum61(n) + mult64(n);	
		sum63(n) = sum62(n) + mult65(n);
		sum64(n) = sum63(n) + mult66(n);	
		sum65(n) = sum64(n) + mult67(n);
		sum66(n) = sum65(n) + mult68(n);	
		sum67(n) = sum66(n) + mult69(n);
		sum68(n) = sum67(n) + mult70(n);	
		sum69(n) = sum68(n) + mult71(n);
		y1(n) = sum69(n) + mult72(n);
	
		
        y(n) = acc;
    end

end
