% add path
addpath(genpath('gloptipoly3'))
addpath('SeDuMi_1_3')

% define constant matrix A
% The nearest normal matrix of A is 
% N = [1.1449+1i*0.8324, -2.0841-1i*0.9958;
% -1.0695 -1i*2.0473,-0.1948- 1i*0.4603]
A = [0.7616+1i*1.2296, -1.4740-1i*0.4577;
    -1.6290-1i*2.6378, 0.1885-1i*0.8575];
A_real = real(A);
A_image= imag(A);
n = size(A,1);


% declare variables
mpol('X_real',n,n)
mpol('X_image',n,n)
LMI = 2; % order of SDP relaxation

P = msdp(min(trace(X_real'*X_real)+trace(X_image'*X_image)),...
    (A_real+X_real)'*(A_real+X_real)+(A_image+X_image)'*(A_image+X_image)-(A_real+X_real)*(A_real+X_real)'-(A_image+X_image)*(A_image+X_image)'==0, ...
    (A_image+X_image)*(A_real+X_real)'-(A_real+X_real)*(A_image+X_image)'+(A_image+X_image)'*(A_real+X_real)-(A_real+X_real)'*(A_image+X_image)==0, ...
    LMI);
[status,obj] = msol(P);

A + double(X_real) +1i*double(X_image)