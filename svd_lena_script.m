% Image Compression with Singular Value Decomposition (SVD).
%   This script uses the SVD for Image Compression, analyses the algorithm
%   (also with Information Theory) and visualizes the results.

close all; clear; clc;


%% Compression

% Original image matrix
Lena_org = imread('Lena.bmp'); % in uint8
Lena = double(Lena_org); % in double

% Call compressing function
compr = 0.01; % change compr to change quality
Lena_red = uint8(svd_compress(Lena,compr));

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
    error('Incorrect input arg: compr must satisfy 0 <= compr <= number of Singular Values');
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
fprintf('Out of %d SVs, only %d SVs saved.\n',r,r_red);
fprintf('Maximum number of SVs for compression: %d SVs.\n',r_max);
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
entropy_red = entropy(Lena_red);
entropy_err = entropy(errorImage);
fprintf('Entropies: Original %d, compressed %d, error %d.\n',entropy_org,entropy_red,entropy_err);

% 2D Histogram: Joint PDF
[jointPDF,~,~] = histcounts2(Lena_org,Lena_red,[m n],'Normalization','probability');

% Joint Entropy
p_logp_nan = jointPDF.*log2(jointPDF);
p_logp = p_logp_nan(isfinite(p_logp_nan));
joint_entropy = -sum(p_logp);
fprintf('Joint entropy: %d.\n',joint_entropy);

% Mutual Information
mi = entropy_org + entropy_red - joint_entropy;
fprintf('Mutual information: %d.\n',mi);


%% Relationship between selcted SVs and made error

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


%% Figures

% Figure 1
figure('units','normalized','outerposition',[0 0 1 1]);

% Original image
subplot(2,3,1)
imshow(uint8(Lena))
title('Original image')

% Histogram of original image
subplot(2,3,4)
imhist(Lena_org);
title('Histogram of original image')

% Compressed image
subplot(2,3,2)
imshow(uint8(Lena_red))
title('Compressed image')

% Histogram of compressed image
subplot(2,3,5)
imhist(Lena_red);
title('Histogram of compressed image')

% Error image
subplot(2,3,3)
imshow(uint8(errorImage))
title('Error image')

% Histogram of error image
subplot(2,3,6)
imhist(errorImage);
title('Histogram of error image')

% Figure 2
figure('units','normalized','outerposition',[0 0 1 1]);

% Compression error over saved SVs
plot(numSVals, displayedError)
xlabel('Number of saved Singular Values')
ylabel('Compression error')
title('Compression error over saved SVs')

% Figure 3
figure('units','normalized','outerposition',[0 0 1 1]);

% 2D Histogram: Joint PDF
histogram2(Lena_org,Lena_red,[m n],'Normalization','probability','FaceColor','flat')
colorbar
