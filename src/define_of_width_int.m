function width = define_of_width_int(a)

width = zeros(length(a),1);

    for i = 1:length(a)
        for j = 1:64
            if (abs(a(i)) < 2^j)
                if (a(i)) < 0
                    width(i) = j+1;
                    break;
                else
                    width(i) = j;
                    break;
                end
            elseif (abs(a(i)) == 2^j)
                width(i) = j+1;
            end
         end
      end

end