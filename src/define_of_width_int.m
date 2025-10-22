function width = define_of_width_int(a)

width = int32(zeros(length(a),1));


    % for i = 1:length(a)
    %     for j = 0:63
    %         if (a(i) < 2^j | a(i) == 2^j)
    %             width(i) = j+1;
    %             break;
    %         end
    %      end
    % end

    for i = 1:length(a)  
        c = a(i) / 2;
        for k = 1:64
            if c > 1
                width(i) = width(i) + 1;
                c = c / 2;
            elseif c == 1
                width(i) = width(i) + 1;
                break;
            end
        end
    end

    width = width + int32(1);
end