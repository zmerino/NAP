function flag = isbeta(str_vec)
    str = ["Beta-a0p5-b1p5","Beta-a2-b0p5","Beta-a0p5-b0p5"];
    flag = zeros(1,length(str_vec));
    for a = 1:size(str_vec,2)
        for b = 1:size(str,2)
            if strcmp(str_vec(a),str(b))
                flag(a) = 1;
                break
            else
                flag(a) = 0;
            end
        end
    end
end