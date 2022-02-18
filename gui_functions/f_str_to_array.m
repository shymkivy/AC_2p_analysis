function array_out = f_str_to_array(str_in)
% inputs
% 'x y z'
% 'x,y,z'
% 'x:y:z'

array1 =  cell(1, numel(str_in));

num_chars = numel(str_in);
start1 = 1;
end1 = 1;
done1 = 0;
while end1 <= num_chars
    
    % if special char
    if end1==start1
        if or(or(str_in(start1:end1) == ' ', str_in(start1:end1) == ','),str_in(start1:end1) == ':')
            array1{start1} = str_in(start1:end1); 
            start1 = end1+1;
            end1 = start1;
            done1 = 1;
        end
    end
    
    % if not end
    if ~done1
        if (end1+1) <= num_chars
            if or(and(str_in(end1+1)>='0',str_in(end1+1)<='9'), str_in(end1+1)=='.')
                end1 = end1 + 1;
                done1 = 1;
            end
        end
    end
     
    if ~done1
        array1{start1} = str2double(str_in(start1:end1));
        start1 = end1+1;
        end1 = start1;
    end
    done1 = 0;
end

remove_el = false(numel(array1),1);
for n_el = 1:numel(array1)
    if isempty(array1{n_el})
        remove_el(n_el) = 1;
    end
    if or(array1{n_el} == ' ', array1{n_el} == ',')
        remove_el(n_el) = 1;
    end
end   
array1(remove_el) = [];

done1 = 0;
created_rep = 0;
while ~done1
    col_loc = strcmpi(array1, ':');
    if sum(col_loc)
        locs1 = find(col_loc);
        if sum(col_loc)>1
            if (locs1(2)-locs1(1)) == 2
                rep_seq = array1{locs1(1)-1}:array1{locs1(1)+1}:array1{locs1(1)+3};
                rep_start = locs1(1)-1;
                rep_end = locs1(1)+3;
                created_rep = 1;
            end
        end
        if ~created_rep
            rep_seq = array1{locs1(1)-1}:array1{locs1(1)+1};
            rep_start = locs1(1)-1;
            rep_end = locs1(1)+1;
            created_rep = 1;
        end
        if created_rep
            array1 = [array1(1:rep_start) rep_seq array1(rep_end:end)];
            created_rep = 0;
        end
    else
        done1 = 1;
    end
end

array_out = cat(2,array1{:});

end