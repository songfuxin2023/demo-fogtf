% Create synthesis data DC1 and Dict A
function [A, Yimg, Xtrue] = createDC1(noise_level)
  %%
    % number of end members
    p = 5;  % fixed for this demo

    %SNR in dB
    SNR = noise_level; 
    % noise bandwidth in pixels of the noise  low pass filter (Gaussian)
    bandwidth = 10000; % 10000 == iid noise
    %bandwidth = 5*pi/224; % colored noise 

    rng(10)

    %% 
    % pure pixels
    x1 = eye(p);

    % mixtures with two materials
    x2 = x1 + circshift(eye(p),[1 0]);

    % mixtures with three materials
    x3 = x2 + circshift(eye(p),[2 0]);

    % mixtures with four  materials
    x4 = x3 + circshift(eye(p),[3 0]);

    % mixtures with four  materials
    x5 = x4 + circshift(eye(p),[4 0]);


    % normalize
    x2 = x2/2;
    x3 = x3/3;
    x4 = x4/4;
    x5 = x5/5;


    % background (random mixture)
    %x6 = dirichlet(ones(p,1),1)';
    x6 = [0.1149 0.0741  0.2004 0.2055, 0.4051]';   % as in the paper

    % build a matrix
    xt = [x1 x2 x3 x4 x5 x6];


    % build image of indices to xt
    imp = zeros(3);
    imp(2,2)=1;


    imind = [imp*1  imp*2 imp* 3 imp*4 imp*5;
        imp*6  imp*7 imp* 8 imp*9 imp*10;
        imp*11  imp*12 imp*13 imp*14 imp*15;
        imp*16  imp*17 imp* 18 imp*19 imp*20;
        imp*21  imp*22 imp* 23 imp*24 imp*25];

    imind = kron(imind,ones(5));

    % set backround index
    imind(imind == 0) = 26;

    % generare frectional abundances for all pixels
    [nl,nc] = size(imind);
    np = nl*nc;     % number of pixels
    for i=1:np
        X(:,i) = xt(:,imind(i));
    end

    %  Size of the image
    nl = size(imind,1);
    nc = size(imind,2);

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


    %% select p endmembers  from A

    % angles (a_1,a_j) \sisizemeq min_angle)
    % supp = 1:p;
    supp = [2 3 4 5 6]; % dont get two Jarosites

    % % Sample endmembers at random
    % supp = randsample(size(A,2), p);

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

    % observed spectral vector
    Y = M*X + noise;
    %

    % create  true X wrt  the library A
    n = size(A,2);
    N = nl*nc;
    XT = zeros(n,N);
    XT(supp,:) = X; 
    Yimg = reshape(Y', nl, nc, L);   
    Xtrue = XT;   
 
 
end