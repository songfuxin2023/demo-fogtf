function [X] = fogtf_solver(M,Y,Hx,Hy,lambda1,lambda2,mu,featureNums)
%%
%   min  0.5*||Y-AX||^2_F + lambda1||W2.*U||_1 + lambda2||Z||_1 lR+(V)
%   s.t. W1X=U,X=V,VH=Z. 

disp({'FoGTF-HU is calculating...'});
% mixing matrix size
[LM,n] = size(M); %  LM  is band numbers ; n  is  endmember numbers
% data set size
[LY,N] = size(Y); % LY is band numbers ; N is pixel numbers ;

%% all initializations

% initializations X0
[UF,SF] = svd(M'*M);
sF = diag(SF);
IF = UF*diag(1./(sF+mu))*UF';
X =IF*M'*Y;

% approximate processing IFv
IFv=( Hx*Hx'+Hy*Hy'+eye(N))\eye(N);
[u,d,v]=svd(IFv);
IFvu=u(:,1:featureNums);
d2=diag(d);
d1=diag(d2(1:featureNums));
v1=v(:,1:featureNums);
IFvd=d1*v1';

% initializations weighted matrix W
W1 = diag(ones(n,1));

% L1norm _ U
reg(1) = 1;
U = W1*X;
D = zeros(n,N);

% Possitivity _ V
reg(2) = 2;
V = X;
E = zeros(n,N);

% Regularization GTF _ U
reg(3) = 3;
Zx= V*Hx;  % horizontal 
Fx = zeros(n,N);

reg(4) = 4;
Zy= V*Hy; % vertical 
Fy = zeros(n,N);


%% calculate objective function 

tol = sqrt(N)*1e-5; 
i=1;
j=1;
res1 = inf;
res2 = inf;
epl = 1e-4;
AL_iters1=10;
AL_iters2=20;

while (j <= AL_iters2) && (sum(abs(res2)) > tol)
    
    % update W
    wu = sqrt(sum( X.^2, 2));
    ws = 1./(wu+ epl);
    W1 = diag(ws);
    W2 = 1./(abs(X)+epl);
    IFx = inv(M'*M + mu.*(eye(n) + W1'*W1)); 
    
while (i <= AL_iters1) && (sum(abs(res1)) > tol)
    
    % update X
    X = IFx*(M'*Y + mu.*(W1'*U-W1'*D+V-E));
    
    % update U
    U = soft(W1*X+D,lambda1/mu.*W2);

    %  update  V
    V =(X + E +(Zx- Fx)*Hx'+(Zy-Fy)*Hy')*IFvu*IFvd; % IFv=IFvu*IFvd
    V = max(V,0);
       
    % update Z (GTF)
    Zx = soft(V*Hx + Fx,lambda2/mu);  % horizontal 
    Zy = soft(V*Hy + Fy,lambda2/mu);  % vertical   

    % update Lagrange multipliers   
    D = D + W1*X - U;
    E = E + X - V;
    Fx = Fx + V*Hx - Zx;
    Fy = Fy + V*Hy - Zy;
       
    % compute residuals
    if mod(i,10) == 1
        st = [];        
        res1(1) = norm(W1*X-U,'fro');
        st = strcat(st,sprintf(' res(%i) = %2.6f',reg(1),res1(1) ));
        res1(2) = norm(X-V,'fro');
        st = strcat(st,sprintf('  res(%i) = %2.6f',reg(2),res1(2) ));    
        res1(3) = norm(V*Hx-Zx,'fro');
        st = strcat(st,sprintf('  res(%i) = %2.6f',reg(3),res1(3) )); 
        res1(4) = norm(V*Hy-Zy,'fro');
        st = strcat(st,sprintf('  res(%i) = %2.6f',reg(4),res1(4) )); 
        fprintf(strcat(sprintf('iter = %i -',j),st ,'\n'));
    end
        
i=i+1;
   
end  

    i=1;
    j=j+1;
    res2(j)=res1(1);
    
end   
end