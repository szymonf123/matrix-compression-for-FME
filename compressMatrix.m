function [v] = compressMatrix(A, epsilon, r)
    
    [numRows,numCols]=size(A);
    
    %if numRows == numCols
    %    n=numRows;
    %else
    %    warning('This is not square matrix')
        %return
    %end
    
    %Make empty node
    v = Node;
    
    %Extract sub-matrices from A
    %A11=A(1:floor(n/2),1:floor(n/2));                     %what if n%2!=0??? Warning: Integer operands are required for colon operator when used as index.
    %A12=A(1:floor(n/2),floor(n/2+1):n); 
    %A21=A(floor(n/2+1):n,1:floor(n/2));
    %A22=A(floor(n/2+1):n,floor(n/2+1):n);
    
    A11=A(1:floor(numRows/2),1:floor(numCols/2));
    A12=A(1:floor(numRows/2),floor(numCols/2+1):numCols); 
    A21=A(floor(numRows/2+1):numRows,1:floor(numCols/2));
    A22=A(floor(numRows/2+1):numRows,floor(numCols/2+1):numCols);
    
    node = handleSubmatrix(A11, epsilon, r);
    append_child(v, node);
    
    node = handleSubmatrix(A12, epsilon, r);
    append_child(v, node);
    
    node = handleSubmatrix(A21, epsilon, r);
    append_child(v, node);
    
    node = handleSubmatrix(A22, epsilon, r);
    append_child(v, node);
end


function [node] = handleSubmatrix(A, epsilon, r)
    if nnz(A) == 0                          %nnz - number of non-zero elements in a matrix
        node=Node;
        node.rank=0;
        [numRows,~]=size(A);
        node.rowsWithZero = numRows;
    else
        [numRows,numCols]=size(A);
        minimum=min(numRows, numCols);
        [U,D,V]=svds(A,minimum);
        eigenvalues = diag(D);
        rank = sum(abs(eigenvalues) > epsilon);

        
        
        if rank == 0
            node=Node;
            node.rank=0;
            node.rowsWithZero = numRows;
        elseif ((rank <= r) && (rank < (numRows/2))) || (rank == 1)
            node=Node;
            node.rank=rank;
            node.eigenvalues = eigenvalues(1:rank);
            %node.V_columns = V(:, 1:rank);
            %node.U_rows = U(1:rank, :);
            node.U_columns = U(:,1:rank);
            V_transposed = V';
            node.V_rows = V_transposed(1:rank,:);   %V transponowane

            for i=1:rank
                node.V_rows(i,:) = node.V_rows(i,:) * node.eigenvalues(i);
            end
        else
            node = compressMatrix(A, epsilon, r);
        end

    end
end