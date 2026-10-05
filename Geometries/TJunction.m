clear all;
%Creating variables
L = 0;
rho = 1;
armRadius = 2.54;
Q = 5;
offsetX = 2.54;
offsetY = 2.54;
phi = linspace(0, 3*pi/4, 10000);

%Regular (CCW) Right Arm
for i = 1:1:length(phi)
    x(i) = armRadius*sin(phi(i))+ offsetX;
    y(i) = -armRadius*cos(phi(i))+ offsetY;
end

%Regular (CCW) Left Arm
for i = 1:1:length(phi)
    xMirror(i) = -armRadius*sin(phi(i))- offsetX;
    yMirror(i) = armRadius*cos(phi(i))- offsetY;
end

for i = (0.3*100):1:(offsetX*100)
    xStraight(i-29) = i/100;
end
yStraight = zeros(1, length(xStraight));

for i = 1:1:630
    xCircle(i) = 0.75*2.54*cos(i/100);
    yCircle(i) = 0.75*2.54*sin(i/100);
end

for i = 1:1:630
    xSmallCircle(i) = 0.15*2.54*cos(i/100);
    ySmallCircle(i) = 0.15*2.54*sin(i/100);
end


%Calculate Angular Momentum from (Flow Rate multiplied by two arms, but divided because flow originally split into two arms)
for i = 1:1:length(x)-1
    ds(i) = sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2);
    r = [x(i) y(i) 0];
    T = [x(i)-x(i+1) y(i)-y(i+1) 0];
    sinVector = cross(r,T)/(norm(r)*norm(T));
    sinTheta(i) = sinVector(3); % REMEMBER: NORM() IS THE PROBLEM CAUSING NO NEGATIVE ANGLES
    
    R(i) = sqrt(x(i)^2+y(i)^2);
    %        𝜌    Q   sin(𝜃)          R     ds
    dL(i) = rho * Q * sinTheta(i) * R(i) * ds(i);
    L = L + dL(i); % Integrating
end

figure('DefaultAxesFontSize',20)

dL(10000) = NaN;
ds(10000) = NaN;
dLStraight = zeros(1, length(xStraight))

cmap = "jet";
R = [linspace(1,1,100) linspace(1, 0, 100)];  %// Red from 1 to 0
B = [linspace(0,1,100) linspace(1, 1, 100)];  %// Blue from 0 to 1
G = [linspace(0,1,100) linspace(1, 0, 100)];   %// Green all zero

colormap( [R(:), G(:), B(:)] );
caxis([-67 67])

hold on
%Outline Arms
plot(x, y, color = [0 0 0], LineWidth=19)
plot(xMirror, yMirror, color = [0 0 0], LineWidth=19)
plot(xStraight, yStraight, color = [0 0 0], LineWidth=19)
plot(-xStraight, -yStraight, color = [0 0 0], LineWidth=19)

patch(xStraight,yStraight,dLStraight,'EdgeColor','interp','Marker','o','MarkerFaceColor','flat', 'MarkerSize', 15);
patch(-xStraight,-yStraight,dLStraight,'EdgeColor','interp','Marker','o','MarkerFaceColor','flat', 'MarkerSize', 15);
patch(x,y,dL./ds,'EdgeColor','interp','Marker','o','MarkerFaceColor','flat', 'MarkerSize', 15);
patch(xMirror,yMirror,dL./ds,'EdgeColor','interp','Marker','o','MarkerFaceColor','flat', 'MarkerSize', 15);

plot(xSmallCircle, ySmallCircle, color = [0.7 0.7 0.7], LineWidth = 5)
plot(xCircle, yCircle, color = [0.7 0.7 0.7], LineWidth = 5, LineStyle="--")
%plot(x, y, color = [0 0 0], LineWidth= 2, LineStyle= ":")
%plot(xStraight, yStraight, color = [0 0 0], LineWidth= 2, LineStyle= ":")

%xlabel("cm")
%ylabel("cm")
%title("Regular (CCW), Re = 600")

box off
axis equal
axis off

%a = colorbar;
%a.Label.String = "Angular Momentum Density (g-cm/s)";


axis equal
fprintf(string(L));
hold off