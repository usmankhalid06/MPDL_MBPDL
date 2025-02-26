function G = gaussian2D(rows, cols, sigma)
    [X, Y] = meshgrid(linspace(-cols/2, cols/2, cols), linspace(-rows/2, rows/2, rows));
    G = exp(-(X.^2 + Y.^2) / (2 * sigma^2));
    G = G / max(G(:));  % Normalize to have a peak of 1
end
