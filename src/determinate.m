
function sum_deet = determinate(a)

% a = [4 6 -2 4 3; 1 2 -3 1 5; 4 -2 1 0 3; 6 4 4 6 8; 3 2 8 1 2];
m = [1 -1 1 -1 1];
% a = [ 4 6  -2 4 3 ] a11 a12 a13 a14 a15
%     [ 1 2  -3 1 5 ] a21 a22 a23 a24 a25
%     [ 4 -2  1 0 3 ] a31 a32 a33 a34 a35
%     [ 6 4   4 6 8 ] a41 a42 a43 a44 a45
%     [ 3 2   8 1 2 ] a51 a52 a53 a54 a55

s.a1 = filloutliers(a(2:5,[2 3 4 5]), 'linear','grubbs'); % delete 1 column
s.a2 = filloutliers(a(2:5,[1 3 4 5]), 'linear','grubbs'); % delete 2 column
s.a3 = filloutliers(a(2:5,[1 2 4 5]), 'linear','grubbs'); % delete 3 column
s.a4 = filloutliers(a(2:5,[1 2 3 5]), 'linear','grubbs'); % delete 4 column
s.a5 = filloutliers(a(2:5,[1 2 3 4]), 'linear','grubbs'); % delete 4 column

x={'a1' 'a2' 'a3' 'a4' 'a5'};

for i = 1:5

    aa1 = filloutliers(s.(x{i})(2:4,[2 3 4]), 'linear','grubbs'); % delete 1 column
    aa2 = filloutliers(s.(x{i})(2:4,[1 3 4]), 'linear','grubbs'); % delete 2 column
    aa3 = filloutliers(s.(x{i})(2:4,[1 2 4]), 'linear','grubbs'); % delete 3 column
    aa4 = filloutliers(s.(x{i})(2:4,[1 2 3]), 'linear','grubbs'); % delete 4 column

    ar1 =  s.(x{i})(1,1) * det_3x3(aa1);
    ar2 =  s.(x{i})(1,2) * det_3x3(aa2);
    ar3 =  s.(x{i})(1,3) * det_3x3(aa3);
    ar4 =  s.(x{i})(1,4) * det_3x3(aa4);

    deet(i) = (m(i)*a(1,i)) * (ar1 - ar2 + ar3 - ar4)
end


sum_deet = 0;
for i = 1:5
    sum_deet = sum_deet + deet(i);
end



end

function [aa] = det_3x3(a2);

    a31 = a2(2:3,2:3);
    a32 = a2(2:3,1:2:end);
    a33 = a2(2:3,1:2);

    a331 = a2(1,1) *(a31(1,1)*a31(2,2) - a31(2,1)*a31(1,2));
    a332 = a2(1,2) *(a32(1,1)*a32(2,2) - a32(2,1)*a32(1,2));
    a333 = a2(1,3) *(a33(1,1)*a33(2,2) - a33(2,1)*a33(1,2));
    aa = a331 - a332 + a333;

end