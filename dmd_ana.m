%% preprocess
clear, clc, clf, close all;

% load regular_plastic.mat;
% load composite.mat;
% load carrier_d1.mat;
load carrier_d2.mat;

data = Fri_coeff;
figure(1);
subplot(2, 3, 1)
plot(Time, Torque)
title('torque')

subplot(2, 3, 2)
plot(Time, Load)
title('load')

subplot(2, 3, 3)
plot(Time, Fri_coeff)
title('fri coefficient')

subplot(2, 3, 4)
plot(Time, Temp)
title('temperature')

subplot(2, 3, 5)
plot(Time, Rotation)
title('rotationspeed')

subplot(2, 3, 6)
plot(Time, Revolution)
title('revolution')

h = 20;
step = 1;

CU = [];
for i = 1:floor((length(Time)-h)/step) - 1
    temp = data(step * i:step * i + h - 1, :);
    CU = [CU, temp(:)];
end

[nx, ny] = size(temp);
X = CU(:,1:end-1);
X2 = CU(:,2:end);
[U,S,V] = svd(X,'econ');

%%  Compute DMD (Phi are eigenvectors)
r = 30;
r = min(r, h);
U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);
Atilde = U'*X2*V*inv(S);
[W,eigs] = eig(Atilde);
Phi = X2*V*inv(S)*W;

%% Plot DMD modes
figure(2);
for i = 1:1:20
    subplot(4, 5, (i-1)/1 + 1);
    plot(abs(Phi(:,i)));
    title(['mode' num2str(i)]);
end
    
%%  Plot DMD spectrum
figure(3)
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'k--') % plot unit circle
hold on, grid on
scatter(real(diag(eigs)),imag(diag(eigs)),'ok')
axis([-1.1 1.1 -1.1 1.1]);

