X1 = 1.5;
X2 = 2;
u1 = @(X3) -X1^0.8 * ((1/(1+X1^0.8 +X3^0.8))+(1/(1+X1^0.8+X2^0.8)));
u2 = @(X3) -X2^0.8 * ((1/(1+X2^0.8 +X3^0.8))+(1/(1+X1^0.8+X2^0.8)));

grid =0:0.05:10;
u1_list = grid;
u2_list = grid;
p1 = grid;
p2 = grid;
for i=1:numel(grid)
    X3 = grid(i);
    u1_list(i) = u1(X3);
    u2_list(i) = u2(X3);
    p1(i) = exp(u1_list(i))/(exp(u1_list(i))+exp(u2_list(i)));
    p2(i) = exp(u2_list(i))/(exp(u1_list(i))+exp(u2_list(i)));
end
%%
plot(grid,u1_list,grid,u2_list)
plot(grid,[p1;p2])

%%
K = 2;
J = 3;
ChoiceSet = @(X31) [1 1;10 10; X31 0.5]; % 3 decoys 1 when ptice below 1
Theta = [1 0 1 1];
attrSign = [-1 1];
% debug
Theta = [1, 0, 5, 7];
% X_orig = Xs{obs}([1 3 2],:); %
X_orig = [17 0.9; 27 1.3; 18 0.6]
X_orig = [43 0.935; 27 1.3; 0 0.6]
attrSign = [1 1];
X_orig(3,1) = 0;
ChoiceSet = @(xdec) X_orig + [0 0; 0 0; xdec 0];

grid =0:0.1:45;
p1 = grid;
p2 = grid;
p3 = grid;
p1_d = grid;
p2_d = grid;
v_list = nan(3,numel(grid));
v_list_d = nan(2,numel(grid));
norm_coefs1 = nan(3,numel(grid));
norm_coefs1_d = nan(3,numel(grid));
for i=1:numel(grid)
    X = ChoiceSet(grid(i));
    X_Dual = X(1:2,:);
     %True params
    alpha = Theta(1);
    sigma = Theta(2);
    Beta = (attrSign .* Theta(3:3+K-1))';
    %utility computation
    u_x = X.^alpha;
    v = zeros(J,1);
    unnorm_u =  Beta .* u_x';
    u_x_d = X_Dual.^alpha;
    v_d = zeros(2,1);
    unnorm_u_d =  Beta .* u_x_d';
    for j=1:J
        u_y = u_x;
        u_y(j,:)=[];
        pairwise_sums = sigma +  u_x(j,:) + u_y
        norm_coefs = sum(1 ./ pairwise_sums,1);%./(J-1);
%         norm_coefs = norm_coefs.^attrSign;
%         norm_coefs(attrSign == -1) = J - norm_coefs(attrSign == -1);
        v(j) = norm_coefs * unnorm_u(:,j);
        if j<=2
            u_y_d = u_x_d;
            u_y_d(j,:)=[];
            norm_coefs_d = sum(1 ./ (sigma + u_x_d(j,:) +  u_y_d),1);
            v_d(j) = norm_coefs_d * unnorm_u_d(:,j);
        end
    end
%     v = v - max(v); %avoid overflow
    v_list(:,i) = v;
    sum_exp_v = sum(exp(v));
    probachoice = exp(v)./sum_exp_v;
    v_list_d(:,i) = v_d;
    sum_exp_v_d = sum(exp(v_d));
    probachoice_d = exp(v_d)./sum_exp_v_d;
    p1(i) = probachoice(1);
    p2(i) = probachoice(2);
    p3(i) = probachoice(3);
    p1_d(i) = probachoice_d(1);
    p2_d(i) = probachoice_d(2);
end
f = figure('Name','Positive','Position',[500 50 700 900]);
subplot(2,1,1);
norm_proba = [p1;p2] ./ (p1+p2);
plot(grid,[p1;p2] ./ (p1+p2));
hold on;
plot(grid,[p1_d;p2_d]);
subplot(2,1,2);
plot(grid,v_list');

%%
subplot(3,1,1);
plot(grid,[p1_d;p2_d]);
title('p1 and p2 (binary case)')
subplot(3,1,2);
plot(grid,[p1;p2] ./ (p1+p2));
title('p1 and p2 (with decoy, x_d \in (0,45))')
subplot(3,1,3);
plot(grid,[p1;p2] ./ (p1+p2));
hold on;
plot(grid,[p1_d;p2_d]);
title('Superimposition')

