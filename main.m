B1 = [1,1,1,1,0,0,0,0;1,2,3,4,0,0,0,0;0,1,2,3,0,0,0,0;2,0,5,6,0,0,0,0;1,1,1,0,0,0,0,0;0,2,0,4,0,2,0,0;0,1,2,0,0,0,0,1;2,0,0,6,0,1,0,0]
[numRows,numCols]=size(B1);
z=rand(numRows,1);
firstCost = nnz(B1)*2
firstResult = B1*z
epsilon=0.01;
r=6;
v = compressMatrix(B1, epsilon, r);
[secondResult, secondCost] = MultiplyMatrixByVector(v, z);
difference = secondResult - firstResult;
generateBitmap(640,v);
normalised_difference=norm(difference)
cost_difference=firstCost-secondCost

