
function C = reshape_groups(groups, data)
    [row, col] = size(data);
    
    group = zeros(col, 1);
    C = [];
    
    for i = 1:row
        new = reshape(data(i,:),[col,1]);
    
        for j = 1:col

            disp(group(j))
            disp(groups(i));
            group(j) = groups(i);
        end
        disp(group)

        C = [C; horzcat(group, new)];
    end
