function [Dh,Dv]= fogdo(I,lambda,sigma,sharpness,maxIter)

    I = im2double(I);
    x = I;
    sigma_iter = sigma;
    lambda = lambda/2.0;
    dec=2.0;
    for iter = 1:maxIter
        [wx, wy, ~, ~] = computeTextureWeights(x, sigma_iter, sharpness);
        [x, ~, ~,Dh,Dv] = solveLinearEquation_in_rog(I, wx, wy, lambda);
       
        sigma_iter = sigma_iter/dec;
        if sigma_iter < 0.5
            sigma_iter = 0.5;
        end
    end
end



function [retx, rety, fx, fy] = computeTextureWeights(fin, sigma,sharpness)   
   fx = diff(fin,1,2);
   fx = padarray(fx, [0 1 0], 'post');
   fy = diff(fin,1,1);
   fy = padarray(fy, [1 0 0], 'post');

   vareps_s = sharpness;
   vareps = 0.0001;

   wto = max(sum((abs(fx)+abs(fy)),3)/size(fin,3),vareps_s).^(-1);
   
   gfx = lpfilter(fx,sigma);
   wtbx = max(sum(abs(gfx),3)/size(fin,3),vareps).^(-1);
   gfy = lpfilter(fy,sigma);
   wtby = max(sum(abs(gfy),3)/size(fin,3),vareps).^(-1);
 
   retx = wtbx.*wto;
   rety = wtby.*wto;

   retx(:,end) = 0;
   rety(end,:) = 0;
end

function ret = conv2_sep(im, sigma)
  ksize = bitor(round(5*sigma),1);
  g = fspecial('gaussian', [1,ksize], sigma);
 
  ret = conv2(im,g,'same');
  ret = conv2(ret,g','same');
end

function FBImg = lpfilter(FImg, sigma)
    FBImg = FImg;
    for ic = 1:size(FBImg,3)
        FBImg(:,:,ic) = conv2_sep(FImg(:,:,ic), sigma);
    end
end

function [out, A, L,Dh,Dv] = solveLinearEquation_in_rog(in, wx, wy, lambda)  
    [h,w,c] = size(in);
    n = h*w;

    wx = wx(:);
    wy = wy(:);
       
    v = ones(n,1);
    wy1 = padarray(wy, 1, 'pre'); wy1 = wy1(1:end-1);
    wx1 = padarray(wx, h, 'pre'); wx1 = wx1(1:end-h);
    Dy = spdiags([-v, v], [-1, 0], n, n); 
    Dx = spdiags([-v, v], [-h, 0], n, n); 
    Wx = spdiags(wx1, 0, n, n); 
    Wy = spdiags(wy1, 0, n, n); 
    L2= Dy'*Wy*Dy+Dx'*Wx*Dx;
    L=sparse(L2); 
    Dx1 = sqrt(Wx)*Dx; 
    Dh = Dx1'; 
    Dy1 = sqrt(Wy)*Dy;
    Dv = Dy1'; 
    
    % solving AS = X
    A = speye(n) + lambda*L; 
    if exist('ichol','builtin')
        F = ichol(A,struct('michol','on'));    
        out = in;
        for i=1:c
            tin = in(:,:,i);
            [tout, ~] = pcg(A, tin(:),0.1,100, F, F');  % AS=X,  tout: S, tin: X.
            out(:,:,i) = reshape(tout, h, w);
        end    
        
    else
        out = in;
        for i=1:c
            tin = in(:,:,i);
            tout = A\tin(:);  % AS=X, tin: X, tout: S
            out(:,:,i) = reshape(tout, h, w);
        end    
    end
end
