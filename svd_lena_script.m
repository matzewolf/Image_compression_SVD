% Image Compression with Singular Value Decomposition (SVD).
%   This script uses the SVD for Image Compression, analyses the algorithm
%   (also with Information Theory) and visualizes the results.

close all; clear; clc;
tic;


%% Compression

% Original image matrix
Lena_org = imread('Lena.bmp'); % in uint8
Lena = double(Lena_org); % in double

% Call compressing function (and measure performance)
compr = 0.01; % change compr to change quality
tic;
Lena_red = uint8(svd_compress(Lena,compr));
func_time = toc;
fprintf('Execution time of svd_compress: %d seconds.\n',func_time);

% Save compressed image
imwrite(Lena_red,'ReducedLena.bmp');


%% Analysis of the algorithm

% SVD on the image
[U,S,V] = svd(Lena);

% Extract Singular Values (SVs)
singvals = diag(S);

% Determine SVs to be saved
if compr >= 0 && compr < 1
    % only SVs bigger than compr times biggest SV
    indices = find(singvals >= compr * singvals(1));
elseif compr >= 1 && compr <= length(singvals)
    % only the biggest compr SVs
    indices = 1:compr;
else
    % return error
    error(...
    'Incorrect input arg: compr must satisfy 0 <= compr <= number of Singular Values');
end

% Size of the image
m = size(Lena,1);
n = size(Lena,2);
storage = m*n;
fprintf('Size of image: %d px by %d px, i.e. uses %d px of storage.\n',m,n,storage);

% SVs and reduced storage
r = min([m,n]); % original number of SVs
r_red = length(indices); % to be saved number of SVs
r_max = floor(m*n/(m+n+1)); % maximum to be saved number of SVs for compression
storage_red = m*r_red + n*r_red + r_red;
if compr >= 0 && compr < 1
    % only SVs bigger than compr times biggest SV
    fprintf('The smallest SV chosen to be smaller than %d of the biggest SV.\n',compr);
elseif compr >= 1 && compr <= length(singvals)
    % only the biggest compr SVs
else
    % return error
    fprintf('There was some error before. Analysis cannot continue.\n')
end
fprintf('Out of %d SVs, only %d SVs saved ',r,r_red);
fprintf('(Maximum number of SVs for compression: %d SVs).\n',r_max);
fprintf('Reduced storage: %d px.\n',storage_red);

% Determine made error
error = 1 - sum(singvals(indices))/sum(singvals);
fprintf('Made error: %d.\n',error);
errorImage = Lena_org - Lena_red;

% Histograms
[dist_org, bin_org] = imhist(Lena_org);
[dist_red, bin_red] = imhist(Lena_red);
[dist_err, bin_err] = imhist(errorImage);

% Entropy
entropy_org = entropy(Lena_org);
fprintf('Entropy of original image: %d bit.\n',entropy_org);
entropy_red = entropy(Lena_red);
fprintf('Entropy of compressed image: %d bit.\n',entropy_red);
entropy_err = entropy(errorImage);
fprintf('Entropy of error image: %d bit.\n',entropy_err);

% 2D Histogram: Joint PDF
[jointPDF,~,~] = histcounts2(Lena_org,Lena_red,[m n],'Normalization','probability');

% Joint Entropy
p_logp_nan = jointPDF.*log2(jointPDF);
p_logp = p_logp_nan(isfinite(p_logp_nan));
joint_entropy = -sum(p_logp);
fprintf('Joint entropy: %d bit.\n',joint_entropy);

% Mutual Information
mi = entropy_org + entropy_red - joint_entropy;
fprintf('Mutual information: %d bit.\n',mi);


%% Relationship between selcted SVs and ...

numSVals = 1:10:r;

% ...used storage
storageSV = m*numSVals + n*numSVals + numSVals;

% ...made error and entropies (compressed and error)
displayedError = zeros(size(numSVals));
entropySV = zeros(4,length(numSVals));
    % 1st row entropy of compressed image, 2nd row entropy of error image
    % 3rd row joint entropy, 4th row mutual information
j = 1; % position in the display vectors
for i = numSVals
    % store S in a temporary matrix
    S_loop = S;
    % truncate S
    S_loop(i+1:end,:) = 0;
    S_loop(:,i+1:end) = 0;
    % construct Image using truncated S
    Lena_red_loop = uint8(U*S_loop*V');
    % construct error image
    Lena_err_loop = Lena_org - Lena_red_loop;
    % compute error
    error_loop = 1 - sum(diag(S_loop))/sum(diag(S));
    % add error to display vector
    displayedError(j) = error_loop;
    % compute entropy of compressed image and add to row 1 of display matrix
    entropySV(1,j) = entropy(Lena_red_loop);
    % compute entropy of error image and add to row 2 of display matrix
    entropySV(2,j) = entropy(Lena_err_loop);
    % compute joint entropy of original and compresed image
    [jointPDF_loop,~,~] = histcounts2(Lena_org,Lena_red_loop,[m n],...
        'Normalization','probability');
    p_logp_nan_loop = jointPDF_loop.*log2(jointPDF_loop);
    p_logp_loop = p_logp_nan_loop(isfinite(p_logp_nan_loop));
    entropySV(3,j) = -sum(p_logp_loop);
    % compute mutual information of original and compressed image
    entropySV(4,j) = entropy_org + entropySV(1,j) - entropySV(3,j);
    % update position
    j = j + 1;
end


%% Figure 1

figure('Name','Images and Histograms','units','normalized','outerposition',[0 0 1 1])

% Original image
subplot(2,3,1)
imshow(uint8(Lena))
title('Original image')

% Histogram of original image
subplot(2,3,4)
imhist(Lena_org)
title('Histogram of original image')

% Compressed image
subplot(2,3,2)
imshow(uint8(Lena_red))
title('Compressed image')

% Histogram of compressed image
subplot(2,3,5)
imhist(Lena_red)
title('Histogram of compressed image')

% Error image
subplot(2,3,3)
imshow(uint8(errorImage))
title('Error image')

% Histogram of error image
subplot(2,3,6)
imhist(errorImage)
title('Histogram of error image')


%% Figure 2

figure('Name','Joint Histogram','units','normalized','outerposition',[0 0 1 1])

% 2D Histogram: Joint PDF
histogram2(Lena_org,Lena_red,[m n],'Normalization','probability','FaceColor','flat')
colorbar
title('Joint Histogram')
zlabel('Joint Probability')


%% Figure 3

figure('Name','Properties over selected Singular Values',...
    'units','normalized','outerposition',[0 0 1 1])

% Used storage over saved SVs
subplot(2,2,1)
plot(numSVals, storage.*ones(size(numSVals))) % original storage (horizontal)
hold on
plot(numSVals, storageSV)
legend('Original storage', 'Storage of SVD','Location','northwest')
xlabel('Number of saved Singular Values')
ylabel('Used storage [px]')
title('Used storage over saved SVs')

% Compression error over saved SVs
subplot(2,2,3)
plot(numSVals, displayedError)
xlabel('Number of saved Singular Values')
ylabel('Compression error [-]')
title('Compression error over saved SVs')

% Entropies over saved SVs
subplot(2,2,[2,4])
plot(numSVals, entropy_org.*ones(size(numSVals))) % original entropy (horizontal)
hold on
plot(numSVals, entropySV)
legend('Original entropy', 'Compression entropy', 'Error entropy',...
    'Joint entropy','Mutual information','Location','southoutside')
xlabel('Number of saved Singular Values')
ylabel('Entropies [bit]')
title('Entropies over saved SVs')


%% Execution time

execution_time = toc;
fprintf('Total execution time of svd_lena_script: %d seconds.\n',execution_time);
