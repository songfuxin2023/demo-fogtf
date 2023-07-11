function demo
    clear; clc;
    
    addpath('dataset');
    addpath('functions');
   
    %% ========= Parameter setting =========

    noise_level = 50; % SNR = 10 20 30 40 and 50
    Dataset = 'DC2';  % DC1 and DC2 
    
    switch Dataset
        case 'DC1'
            [A, Yimg, Xtrue] = createDC1(noise_level); 
            
            paramsRTV.sigma = 1.2;
            paramsRTV.sharpness = 0.001;
            paramsRTV.maxIter = 5;
            paramsRTV.fuseBins = 5;
            featureNums = 70; 
            
            switch noise_level
                case 10
                    paramsRTV.lambda = 0.004;               
                    lambda1 = 2e-3;
                    lambda2 = 1e-3;
                    mu = 0.2;      
                case 20
                    paramsRTV.lambda = 0.001;                    
                    lambda1 = 2e-4;  
                    lambda2 = 3e-4;
                    mu = 0.2;   
                case 30
                    paramsRTV.lambda = 0.001;                  
                    lambda1 = 3e-5;  
                    lambda2 = 1e-4;
                    mu = 0.05;   
                case 40 
                    paramsRTV.lambda = 0.0005;                  
                    lambda1 = 1e-4;  
                    lambda2 = 1e-5;
                    mu = 0.05;   
                case 50
                    paramsRTV.lambda = 0.0001;                  
                    lambda1 = 1e-4;  
                    lambda2 = 1e-6;
                    mu = 0.05;   
            end
        case 'DC2'
            [A, Yimg, Xtrue] = createDC2(noise_level); 

            paramsRTV.sigma = 0.5;
            paramsRTV.sharpness = 0.02;
            paramsRTV.maxIter = 5;
            paramsRTV.fuseBins = 5;
            featureNums = 1500; 
            
            switch noise_level
                case 10                  
                    paramsRTV.lambda = 1e-3;                    
                    lambda1 = 4e-4;  
                    lambda2 = 1e-3;
                    mu = 0.3;
                case 20
                    paramsRTV.lambda = 5e-4;                   
                    lambda1 = 5e-4;
                    lambda2 = 2e-4;
                    mu = 0.07;    
                case 30
                    paramsRTV.lambda = 1e-4;                     
                    lambda1 = 1e-3;
                    lambda2 = 5e-5;
                    mu = 0.006; 
                case 40
                    paramsRTV.lambda = 1e-4;                     
                    lambda1 = 1e-3;
                    lambda2 = 2e-5;
                    mu = 0.001; 
                case 50
                    paramsRTV.lambda = 1e-4;                     
                    lambda1 = 1e-3;
                    lambda2 = 1e-5;
                    mu = 0.0001; 
            end                               
    end

%% ===========  FoGTF method ============ %%

% band, endmember and pixel numbers
n_endmembers=size(A,2);
[n_row,n_col,n_bands]=size(Yimg);
XtruethImg=reshape(Xtrue', n_row,n_col, n_endmembers); 
Y_input =  reshape(Yimg,n_row*n_col,n_bands)';

% calculate firt-order graph difference operator
Ya = average_fusion(Yimg, paramsRTV.fuseBins);     
[no_lines,no_rows,no_bands] = size(Ya);
fimg=reshape(Ya,[no_lines*no_rows no_bands]);
[fimg] = scale_new(fimg);
fimg=reshape(fimg,[no_lines no_rows no_bands]);     
[Dh,Dv]= fogdo(fimg,paramsRTV.lambda,paramsRTV.sigma, paramsRTV.sharpness, paramsRTV.maxIter);

% solve the objective function   
[Xout] = fogtf_solver(A, Y_input,Dh,Dv,lambda1,lambda2,mu,featureNums);
X_img = reshape(Xout', n_row, n_col, n_endmembers);  

% SRE value
SREv = 20*log10(norm(Xtrue,'fro')/norm(Xout - Xtrue,'fro'));
SRE.FoGTF_HU = SREv; 

% RMSE values 
RMSE_X= rmse(Xtrue, Xout);
Rmse.FoGTF_HU = RMSE_X;

disp('SRE:')
disp(SRE);  
disp('RMSE')
disp(Rmse);
    
% generate abundance map 
switch Dataset
    case 'DC1'
       str = sprintf('SRE: %.3f', SREv);
       figure('Name', str);
       subplot(251)
       imagesc(XtruethImg(:,:, 2),[0 1]);
       subplot(252)
       imagesc(XtruethImg(:,:, 3),[0 1]);
       subplot(253)
       imagesc(XtruethImg(:,:, 4),[0 1]);
       subplot(254)
       imagesc(XtruethImg(:,:, 5),[0 1]);
       subplot(255)
       imagesc(XtruethImg(:,:, 6),[0 1]);        
       subplot(256)
       imagesc(X_img(:,:, 2),[0 1]);
       subplot(257)
       imagesc(X_img(:,:, 3),[0 1]);
       subplot(258)
       imagesc(X_img(:,:, 4),[0 1]);
       subplot(259)
       imagesc(X_img(:,:, 5),[0 1]);
       subplot(2,5,10)
       imagesc(X_img(:,:, 6),[0 1]);
       colormap jet;  
    case 'DC2'
       str = sprintf('SRE: %.3f', SREv);
       figure('Name', str);
       subplot(2,9,1)
       imagesc(XtruethImg(:,:, 2),[0 1]);
       subplot(2,9,2)
       imagesc(XtruethImg(:,:, 3),[0 1]);
       subplot(2,9,3)
       imagesc(XtruethImg(:,:, 4),[0 1]);
       subplot(2,9,4)
       imagesc(XtruethImg(:,:, 5),[0 1]);
       subplot(2,9,5)
       imagesc(XtruethImg(:,:, 6),[0 1]);
       subplot(2,9,6)
       imagesc(XtruethImg(:,:, 7),[0 1]);
       subplot(2,9,7)
       imagesc(XtruethImg(:,:, 8),[0 1]);
       subplot(2,9,8)
       imagesc(XtruethImg(:,:, 9),[0 1]);  
       subplot(2,9,9)
       imagesc(XtruethImg(:,:, 10),[0 1]);       
       subplot(2,9,10)
       imagesc(X_img(:,:, 2),[0 1]);
       subplot(2,9,11)
       imagesc(X_img(:,:, 3),[0 1]);
       subplot(2,9,12)
       imagesc(X_img(:,:, 4),[0 1]);
       subplot(2,9,13)
       imagesc(X_img(:,:, 5),[0 1]);
       subplot(2,9,14)
       imagesc(X_img(:,:, 6),[0 1]);
       subplot(2,9,15)
       imagesc(X_img(:,:, 7),[0 1]);
       subplot(2,9,16)
       imagesc(X_img(:,:, 8),[0 1]);
       subplot(2,9,17)
       imagesc(X_img(:,:, 9),[0 1]);    
       subplot(2,9,18)
       imagesc(X_img(:,:, 10),[0 1]);
       colormap jet; 
end                 

    rmpath('dataset');
    rmpath('functions');

end