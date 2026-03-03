riga = 1;
nxx  = 5;
pxx  = 2;
rxx  = 0;
ex = 2;
ey = 1;

% Generowanie macierzy masy
A = masslowrankForRange(riga, nxx, pxx, rxx, ex, ey);

% Podstawowe informacje
disp('Rozmiar macierzy:');
disp(size(A));

disp('Liczba niezerowych elementów:');
disp(nnz(A));

disp('Gęstość macierzy (% niezerowych):');
disp(100 * nnz(A) / numel(A));

% Wizualizacja struktury
figure;
spy(A);
title('Struktura macierzy masy');

% Test mnożenia przez wektor
z = rand(size(A, 1), 1);
b = A * z;

disp('Norma wyniku A*z:');
disp(norm(b));
