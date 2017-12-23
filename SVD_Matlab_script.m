close all
clear all
clc

%% 1. Image of Lena
Lena = double(imread('lena.bmp'));
m = size(Lena,1)
n = size(Lena,2)
storage = m*n

%% 2. SVD on Lena
[U,S,V] = svd(Lena);

%% 3. Extract singular values
singvals = diag(S);
r = length(singvals)

%% 4. Determine truncation of the U, S, V matrices
c = 0.01 % change c to change quality
indices = find(singvals >= c * singvals(1));
r_red = length(indices)
r_max = m*n/(m+n+1)
storage_red = m*r_red + n*r_red + r_red

%% 5. Truncate U, S, V matrices
U_red = U(:,indices);
S_red = S(indices,indices);
V_red = V(:,indices);

%% 6. Low-rank approximation of Lena
Lena_red = U_red * S_red * V_red';

%% 7. Save and print result
imwrite(uint8(Lena_red),'Reduced Lena.bmp');
subplot(1,2,1)
imshow(uint8(Lena))
subplot(1,2,2)
imshow(uint8(Lena_red))

%% 8. Error
error = 1 - sum(singvals(indices))/sum(singvals)
errorImage = Lena - Lena_red;
figure
imshow(uint8(errorImage))


%% REMARK: Relationship between selcted Singular Values and made error:

numSVals = 1:10:r;
displayedError = [];

for i = numSVals
    % store S in a temporary matrix
    S_loop = S;
    % truncate S
    S_loop(i+1:end,:) = 0;
    S_loop(:,i+1:end) = 0;
    % construct Image using truncated S
    Lena_loop = U*S_loop*V';
    % compute error
    error_loop = 1 - sum(diag(S_loop))/sum(diag(S));
    % add error to display vector
    displayedError = [displayedError, error_loop];
end

figure
plot(numSVals, displayedError)
xlabel('Number of saved Singular Values')
ylabel('Compression error')
