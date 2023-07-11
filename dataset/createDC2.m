function [A, Yimg, Xtrue] = createDC2(noise_level)
    %%
    %SNR in dB
    SNR = noise_level; 
    % noise bandwidth in pixels of the noise  low pass filter (Gaussian)
    bandwidth = 10000; % 10000 == iid noise
    %bandwidth = 5*pi/224; % colored noise 

    rng(10);

    %% 
    load spatial2.mat

    %  Size of the images
    nl = size(Xim,1);
    nc = size(Xim,2);
    np = nl*nc;     % number of pixels

    %%
    load USGS_1995_Library.mat
    %  order bands by increasing wavelength
    [dummy index] = sort(datalib(:,1));
    A =  datalib(index,4:end);
    names = names(4:end,:);

    % prune the library 
    % min angle (in degres) between any two signatures 
    % the larger min_angle the easier is the sparse regression problem
    min_angle = 4.44;       
    [A, index] = prune_library2(A,min_angle); % 240  signature 
    names = names(index',:);

    % order  the columns of A by decreasing angles 
    [A, index, angles] = sort_library_by_angle(A);
    names = names(index',:);
    namesStr = char(names);

    % Names of the first 10 ordered materials, with 4.44 deg. prunning:
    % 1 - Jarosite GDS99 K,Sy 200C
    % 2 - Jarosite GDS101 Na,Sy 200
    % 3 - Anorthite HS349.3B 
    % 4 - Calcite WS272 
    % 5 - Alunite GDS83 Na63 
    % 6 - Howlite GDS155
    % 7 - Corrensite CorWa-1
    % 8 - Fassaite HS118.3B  
    % 9 - Adularia GDS57 Orthoclase  
    % 10 - Andradite NMNH113829 

    supp = [2 3 4 5 6 7 8 9 10]; % dont take 2 Jarosites

    M = A(:,supp);
    [L,p] = size(M);  % L = number of bands; p = number of material


    %%
    % set noise standard deviation
    sigma = sqrt(sum(sum((M*X).^2))/np/L/10^(SNR/10));
    % generate Gaussian iid noise
    noise = sigma*randn(L,np);

    % make noise correlated by low pass filtering
    % low pass filter (Gaussian)
    filter_coef = exp(-(0:L-1).^2/2/bandwidth.^2)';
    scale = sqrt(L/sum(filter_coef.^2));
    filter_coef = scale*filter_coef;
    noise = idct(dct(noise).*repmat(filter_coef,1,np));

    %  observed spectral vector
    Y = M*X + noise;

    % create  true X wrt the library A
    n = size(A,2);
    N = nl*nc;
    XT = zeros(n,N);
    XT(supp,:) = X;    
    Yimg = reshape(Y', nl, nc, L);   
    Xtrue = XT;
    
       
end