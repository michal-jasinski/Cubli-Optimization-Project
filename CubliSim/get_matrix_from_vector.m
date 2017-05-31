function matrix = get_matrix_from_vector(vector)
    
    row = length(vector)/3;
    col = 3;
    k = 1;
    matrix = zeros(row,col);
    
    for i = 1 : row
        for j = 1 : col
            matrix(i,j) = vector(k);
            k = k + 1;
        end
    end
end

