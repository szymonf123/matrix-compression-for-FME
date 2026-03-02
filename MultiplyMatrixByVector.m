function [result, cost] = MultiplyMatrixByVector(node, v)
%Description:
%   Function multiply compressed matrix (generated with the usage of function 'compressMatrix') by indicated vector.
%Parameters:
%   node - root node of the compressed matrix
%   v - vector
%Returned variable:
%   result - result of multiplying compressed matrix 'node' by vector 'v'. It means, it is a vector (matrix with one column).

    cost = 0;
    if node.no_of_children == 0
        if node.rank == 0
%             [numRows,~]=size(v);
%             zeroVector = zeros(node.rowsWithZero,numRows);
%             result = zeroVector*v;
              result = zeros(node.rowsWithZero,1); %wektor zer
        else

            result = node.U_columns * (node.V_rows * v);
            cost = cost + nnz(node.V_rows)*2;
            cost = cost + nnz(node.U_columns)*2;
        end
    else
        
        [numRows,~]=size(v);

        v1 = v(1:floor(numRows/2),:);
        v2 = v(floor(numRows/2+1):numRows, :);
        
        %[numRows1,numCols1]=size(v1);
        %[numRows2,numCols2]=size(v2);

        [res1, nested_cost1] = MultiplyMatrixByVector(node.children(1), v1);
        cost = cost + nested_cost1;
        [res2, nested_cost2] = MultiplyMatrixByVector(node.children(2), v2);
        cost = cost + nested_cost2;
        [res3, nested_cost3] = MultiplyMatrixByVector(node.children(3), v1);
        cost = cost + nested_cost3;
        [res4, nested_cost4] = MultiplyMatrixByVector(node.children(4), v2);
        cost = cost + nested_cost4;
        
        if (not(any(res1))) %same zera w res1
           res1res2 = res2;
        elseif (not(any(res2)))
                res1res2 = res1; %same zera w res2
        else
           res1res2 = res1 + res2;
           [numRows12,~]=size(res1res2);
           cost = cost+numRows12;
        end
        
        if (not(any(res3))) %same zera w res3
           res3res4 = res4;
        elseif (not(any(res4)))
                res3res4 = res3; %same zera w res4
        else
           [numRows3,numCol3]=size(res3);
           [numRows4,numCol4]=size(res4);
           res3res4 = res3 + res4;
           [numRows34,~]=size(res3res4);
           cost = cost+numRows34;
        end

        result = [res1res2; res3res4];
    end
end
