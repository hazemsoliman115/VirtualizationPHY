clear all

clc
tot_dist = [];
proj_dist = [];
epsi = 2;
%x=randn(8,1);%+1i*randn(8,1);


prob_vec = zeros(7,1);
z_x = zeros(10,1);
z_y = zeros(10,1);
for m_k = 0:7
    m_k

    proj_dist = zeros(1,100000);
    tot_dist = zeros(1,100000);
    %pd = makedist('Poisson','lambda',m_k);
    %t = truncate(pd,0,7);
    for i=1:100000
    %     H = randn(2,8);%+1i*randn(4,8);
        %H = randn(6,8);%+1i*randn(4,8);
        R = poissrnd(m_k);
        if R >= 8
            continue
        end

        %R = random(t,1,1);
        H = randn(R,8)+1i*randn(R,8);
        x=randn(8,1)+1i*randn(8,1);


        [U,S,V] = svd(H);

        S_t = eye(8,8);
        for s_i=1:R
            S_t(s_i,s_i) = 0;
        end
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
        proj_dist(i) = x'*y;
        tot_dist(i) = x'*x;

    end
    z_x = z_x/100000;
    z_y = z_y/100000;
    prob_vec(m_k+1) = sum(proj_dist > epsi)/100000;
end
%histfit(real(proj_dist),40,'beta')
%dfittool(real(proj_dist).^2)
for m_k = 0:7
    low_bound(m_k+1) = 1/epsi*(8-m_k);
end

plot([0:7], prob_vec,'DisplayName','prob_vec','YDataSource','prob_vec');figure(gcf)
hold on
plot([0:7], low_bound,'r')