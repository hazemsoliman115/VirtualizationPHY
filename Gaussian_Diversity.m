clear all

clc
tot_dist = [];
proj_dist = [];
%x=randn(8,1);%+1i*randn(8,1);
H = randn(6,8)+1i*randn(6,8);

z_x = zeros(10,1);
z_y = zeros(10,1);
for i=1:100000
%     H = randn(2,8);%+1i*randn(4,8);
    %H = randn(6,8);%+1i*randn(4,8);
    x=randn(8,1)+1i*randn(8,1);
    

    [U,S,V] = svd(H);

    S_t = eye(8,8);
    S_t(1,1) = 0;
    S_t(2,2) = 0;
    S_t(3,3) = 0;
    S_t(4,4) = 0;
    S_t(5,5) = 0;
    S_t(6,6) = 0;
    P = V*S_t*V';

    
%     y=P*[x; zeros(size(H,2)-size(x,1),1)];
% 
%     proj_dist = [proj_dist [x; zeros(size(H,2)-size(x,1),1)]'*y];
%     tot_dist = [tot_dist [x; zeros(size(H,2)-size(x,1),1)]'*[x; zeros(size(H,2)-size(x,1),1)]];
    y=P*x;
    y_hat = y;
    %x=x/sqrt(sum(abs(x).^2));
    y=y/sqrt(sum(abs(y).^2));
    for j=1:10
        z_x(j) = z_x(j) + max(abs(x).^2 > j/10);
        z_y(j) = z_y(j) + max(abs(y).^2 > j/10);
    end
    proj_dist = [proj_dist x'*y];
    tot_dist = [tot_dist x'*x];
    
end
z_x = z_x/100000;
z_y = z_y/100000;
%histfit(real(proj_dist),40,'beta')
dfittool(real(proj_dist).^2)