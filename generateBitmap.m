function generateBitmap(N, v)
    %generuje bitmapę N x N
    %v - korzeń drzewa
    %sposób działania: 
    %   - tworzymy macierz NxN przyjmującą wartości 0 lub 1
    %   - 0 - wartości oznaczające podczas wyświetlania czarny piksel, 1 - biały piksel
    
    bitmapMatrix = zeros(N,N);
    
    bitmapMatrix(:) = 1;        %tworzenie białego tła na całej bitmapie
    
    %rysowanie głównej ramki - start
    bitmapMatrix(:, 1) = 0;     %ustawiamy pierwszą kolumnę na 0
    bitmapMatrix(:, end) = 0;   %ustawiamy ostatnią kolumnę na 0
    
    bitmapMatrix(1, :) = 0;     %ustawiamy pierwszy wiersz na 0
    bitmapMatrix(end, :) = 0;   %ustawiamy ostatni wiersz na 0
    %rysowanie głównej ramki - end
    
    
    %dzielenie macierzy na 4 podmacierze - start
    idx_of_first_column = 1;
    idx_of_last_column = N;
    
    idx_of_first_row = 1;
    idx_of_last_row = N;
    
    idx_of_mid_column = N/2;
    idx_of_mid_row = N/2;
    
    bitmapMatrix(idx_of_first_row : idx_of_last_row, idx_of_mid_column) = 0;     %ustawiamy środkową kolumnę na 0
    bitmapMatrix(idx_of_mid_row, idx_of_first_column:idx_of_last_column) = 0;    %ustawiamy środkowy wiersz na 0
    
    %dzielenie macierzy na 4 podmacierze - end
    
    %obsługa każdej z ćwiartek macierzy - start
    bitmapMatrix = handleSubMatrix(bitmapMatrix,v.children(1), idx_of_first_row, idx_of_mid_row, idx_of_first_column, idx_of_mid_column);   %A11
    bitmapMatrix = handleSubMatrix(bitmapMatrix,v.children(2), idx_of_first_row, idx_of_mid_row, idx_of_mid_column, idx_of_last_column);    %A12
    bitmapMatrix = handleSubMatrix(bitmapMatrix,v.children(3), idx_of_mid_row, idx_of_last_row, idx_of_first_column, idx_of_mid_column);    %A21
    bitmapMatrix = handleSubMatrix(bitmapMatrix,v.children(4), idx_of_mid_row, idx_of_last_row, idx_of_mid_column, idx_of_last_column);     %A22
    %obsługa każdej z ćwiartek macierzy - end
    
    imshow(bitmapMatrix);

end

function [matrix] = handleSubMatrix(matrix, node, rows_from, rows_to, columns_from, columns_to)   
    if isempty(node.rank)   %rank is empty - więc dzielimy fragment matrix na kolejne podmacierze i wywołujemy handleSubmatrix dla każdego z dzieci z node
        matrix = divideMatrixOn4Submatrices(matrix, rows_from, rows_to, columns_from, columns_to);

        idx_of_mid_column = floor(columns_from+(columns_to - columns_from)/2);  %ERROR - wyszło 249.5)
        idx_of_mid_row = floor(rows_from +(rows_to - rows_from)/2);

        matrix = handleSubMatrix(matrix, node.children(1),rows_from, idx_of_mid_row, columns_from, idx_of_mid_column);
        matrix = handleSubMatrix(matrix, node.children(2),rows_from, idx_of_mid_row, idx_of_mid_column, columns_to);
        matrix = handleSubMatrix(matrix, node.children(3),idx_of_mid_row, rows_to, columns_from, idx_of_mid_column);
        matrix = handleSubMatrix(matrix, node.children(4),idx_of_mid_row, rows_to, idx_of_mid_column, columns_to);
    elseif node.rank ~= 0
        matrix = drawLines(matrix, node.rank, rows_from, rows_to, columns_from, columns_to); %rysujemy kolumny i wiersze (kreski) o grubości k w danym obszarze
    %else
        %return  - %rank == 0 - więc nic nie rysujemy, zostawiamy pusty kwadrat
    end
end

function [matrix] = divideMatrixOn4Submatrices(matrix, rows_from, rows_to, columns_from, columns_to)
    idx_of_mid_column = floor(columns_from+(columns_to - columns_from)/2);  %ERROR - wyszło 249.5)
    idx_of_mid_row = floor(rows_from +(rows_to - rows_from)/2);
    
    matrix(rows_from : rows_to, idx_of_mid_column) = 0;     %ustawiamy środkową kolumnę na 0
    matrix(idx_of_mid_row, columns_from:columns_to) = 0;    %ustawiamy środkowy wiersz na 0
end

function [matrix] = drawLines(matrix, rank, rows_from, rows_to, columns_from, columns_to)
    %matrix(rows_from:rows_to, columns_from:columns_to) = 0;    %zapełnienie pola czarnym tłem zamiast liniami - wersja uproszczona 
    
    matrix(rows_from+1:rows_from+rank, columns_from+1:columns_to-1) = 0;
    matrix(rows_from+1:rows_to-1, columns_from+1:columns_from+rank) = 0;
%     no_of_pixels = 0;
%     for i = rows_from+1 : rows_to-1
%         if no_of_pixels < rank
%             matrix(i,columns_from+1 : columns_to-1) = 0;
%             no_of_pixels = no_of_pixels+1;
%         elseif no_of_pixels == rank
%             no_of_pixels = 0;
%         end
%     end
% 
%     no_of_pixels = 0;
%     for j = columns_from+1 : columns_to-1
%         if no_of_pixels < rank
%             matrix(rows_from+1 : rows_to-1, j) = 0;
%             no_of_pixels = no_of_pixels+1;
%         elseif no_of_pixels == rank
%             no_of_pixels = 0;
%         end        
%     end
end