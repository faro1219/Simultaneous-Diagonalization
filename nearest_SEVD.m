%% add path
addpath(genpath('gloptipoly3'))
addpath(genpath('YALMIP-master'));

%% define constant matrix A
% Example 1
% randomly generated
n = 2;
A = randn(n,n)+1i*randn(n,n);
B = randn(n,n)+1i*randn(n,n);

A_real = real(A);
A_image= imag(A);
B_real = real(B);
B_image= imag(B);

n = size(A,1);


%% declare variables and solve by SDP relaxation
mpol('X_real',n,n)
mpol('X_image',n,n)
mpol('Y_real',n,n)
mpol('Y_image',n,n)
LMI = 2; % order of SDP relaxation

P = msdp(min(trace(X_real'*X_real)+trace(X_image'*X_image)+trace(Y_real'*Y_real)+trace(Y_image'*Y_image)),...
    [(A_real+X_real)'*(A_real+X_real)+(A_image+X_image)'*(A_image+X_image)-(A_real+X_real)*(A_real+X_real)'-(A_image+X_image)*(A_image+X_image)'==0, ...
    (A_image+X_image)*(A_real+X_real)'-(A_real+X_real)*(A_image+X_image)'+(A_image+X_image)'*(A_real+X_real)-(A_real+X_real)'*(A_image+X_image)==0, ...
    (B_real+Y_real)'*(B_real+Y_real)+(B_image+Y_image)'*(B_image+Y_image)-(B_real+Y_real)*(B_real+Y_real)'-(B_image+Y_image)*(B_image+Y_image)'==0, ...
    (B_image+Y_image)*(B_real+Y_real)'-(B_real+Y_real)*(B_image+Y_image)'+(B_image+Y_image)'*(B_real+Y_real)-(B_real+Y_real)'*(B_image+Y_image)==0, ...
    (A_real+X_real)*(B_real+Y_real)-(A_image+X_image)*(B_image+Y_image)-(B_real+Y_real)*(A_real+X_real)+(B_image+Y_image)*(A_image+X_image)==0, ...
    (A_image+X_image)*(B_real+Y_real)+(A_real+X_real)*(B_image+Y_image)-(B_image+Y_image)*(A_real+X_real)-(B_real+Y_real)*(A_image+X_image)==0
    ], ...
    LMI);
[status,obj] = msol(P);

S = A + double(X_real) + 1i*double(X_image);
T = B + double(Y_real) + 1i*double(Y_image);

S*S'-S'*S
T*T'-T'*T
S*T-T*S