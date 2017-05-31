function vector = get_vector_from_matrix(matrix)
    
    [row, col] = size(matrix);
    vector = zeros(row*col);
    k = 1;
    
    for i = 1 : row
        for j = 1 : col
            vector(k) = matrix(i,j);
            k = k + 1;
        end
    end
end

