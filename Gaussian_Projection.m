clear all

clc
tot_dist = [];
proj_dist = [];
%x=randn(8,1);%+1i*randn(8,1);
x=randn(8,1)+1i*randn(8,1);
x=x/sqrt(sum(abs(x).^2));
for i=1:100000
%     H = randn(2,8);%+1i*randn(4,8);
    %H = randn(6,8);%+1i*randn(4,8);
    H = randn(4,8)+1i*randn(4,8);

    [U,S,V] = svd(H);

    S_t = eye(8,8);
    S_t(1,1) = 0;
    S_t(2,2) = 0;
    S_t(3,3) = 0;
    S_t(4,4) = 0;
    %S_t(5,5) = 0;
    %S_t(6,6) = 0;
    P = V*S_t*V';

    
%     y=P*[x; zeros(size(H,2)-size(x,1),1)];
% 
%     proj_dist = [proj_dist [x; zeros(size(H,2)-size(x,1),1)]'*y];
%     tot_dist = [tot_dist [x; zeros(size(H,2)-size(x,1),1)]'*[x; zeros(size(H,2)-size(x,1),1)]];
    y=P*x;

    proj_dist = [proj_dist x'*y];
    tot_dist = [tot_dist x'*x];
    
end

histfit(real(proj_dist),10,'beta')
dfittool(real(proj_dist))