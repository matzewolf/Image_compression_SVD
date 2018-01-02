function [ A_red ] = svd_compress( A_org, compr )

% svd_compress compresses an input matrix (e.g. an image) using the
% Singular Value Decomposition (SVD).
%   Input args: A_org: Any matrix with double real entries, e.g. an image 
%   file (converted from uint8 to double).
%   compr: Quality of compression. If 0 <= compr < 1, it only keeps
%   Singular Values (SVs) larger than compr times the biggest SV. If 1 <= 
%   compr <= number of SVs, it keeps the biggest compr SVs. Otherwise the 
%   function returns an error.
%   Output args: A_red: Compressed version of A_org in double using the
%   SVD, e.g. an image file (convert from double to uint8).

% SVD on the original matrix
[U,S,V] = svd(A_org);

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

% Truncate U,S,V
U_red = U(:,indices);
S_red = S(indices,indices);
V_red = V(:,indices);

% Calculate compressed matrix
A_red = U_red * S_red * V_red';

end
