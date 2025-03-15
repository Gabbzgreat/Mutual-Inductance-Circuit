function Rlc
h = 0.0001;
T = 0:h:30;
y = zeros(3, length(T));
e = zeros(1, length(T));
M_trace = zeros(1, length(T));
uL1_trace = zeros(1, length(T));
UR2 = zeros(1, length(T));
dy = [0
0
0];
M_prev = fM(0);
R = 10;
R1 = 0.1;
R2 = 10;
P= zeros(1, length(T));
for i = 1:length(T) - 1
%E = 240*sin (T(i)); %checking for different sin
E = 210*sin (2*pi*5*T(i));
%E=1;
% P(i+1)=R1*y(1,i)^2+R2*y(2,i)^2;
% t_rect = mod(T(i),3);
%if t_rect < 1.5
% E = 120;
% else
% E = 0;
% end
e(i) = E;
% y(:,i+1) = y(:,i) + h * f(T(i) , y(:,i),E);
o=dy(2);
[dy,M_prev,uL1] = f(T(i), y(:,i), E, M_prev, dy(2));
M_trace(i+1) =M_prev;
uL1_trace(i+1)=uL1;
UR2(i+1)=y(2,i)*R;
% y(:,i+1) = y(:,i) + h * dy;%Euler
b=dy;
%{%Heuns
% [dy2,M_prev,uL1] = f(T(i+1), y(:,i+1), E, M_trace(i), dy(2));
% y(:,i+1)=y(:,i)+h/2*(b+dy2);
%}
c=(T(i)+T(i+1))/2;
k1 = b;
[k2,M_prev,uL1]= f(c, y(:,i) + h/2*k1,E,M_trace(i),o);
[k3,M_prev,uL1]= f(c, y(:,i) + h/2*k2,E,M_trace(i),o);
[k4,M_prev,uL1]= f(T(i)+h, y(:,i) + h*k3,E,M_trace(i),o);
y(:,i+1) = y(:,i) + h/6*(k1+2*k2+2*k3+k4);
end
close all;
hold on
figure (1);
plot(T, y(3,:))
legend('u_C');
xlabel('t[s]');
ylabel('u_C[V]');
grid;
figure (2);
plot(T, y(1,:), T, y(2,:));
legend('i_1', 'i_2');
xlabel('t[s]');
ylabel('i[A]');
grid;
figure (3)
plot(T,(e))
xlabel('t[s]');
ylabel('Etrace[V]');
legend('Etrace');
grid;
figure (4)
plot(T,M_trace)
xlabel('t[s]');
ylabel('M_trace');
legend('M_trace');
grid;
figure (5)
plot(T,uL1_trace)
xlabel('t[s]');
ylabel('uL1_trace[V]');
legend('uL1_trace');
grid;
end
function [Ud, Md ]= m_data()
Ud = [20, 50, 100, 150, 200, 250, 280, 300]; %table
Md = [0.46, 0.64, 0.78, 0.68, 0.44, 0.23, 0.18, 0.18];
%a = polyfit(Ud, Md,4);
%ux = 0:1:max(Ud)
%mx = polyval(a, ux);
%plot(Ud, Md, 'or', ux, mx, 'b')
global a
a = findnewton(Ud,Md);
end
function y = polynewton(a,xin,x)
y = zeros(length(x), 1);
for i = 1:length(x)
s = 0;
for j=1:length(a)
p = 1;
for k=1:(j-1)
p = p * (x(i) - xin(k));
end
s = s + a(j) * p;
end
y(i) = s;
end
end
function a = findnewton(xin, yin)
n = length(xin);
A = zeros(n,n+1); % this is the matrix representing our calculations
A(:,1) = xin;
A(:,2) = yin;
for j=2:n
for i=1:n-j+1
A(i, j+1) = (A(i+1,j) - A(i,j)) / (xin(i+j-1) - xin(i));
end
end
A;
a = A(1,2:end); % here we slice out part of the first row
end
function [dy, M, uL1] = f(t, y, e, M_prev, dy2_prev)
R1 = 0.1;
R2 = 10;
C = 0.5;
L1 = 3;
L2 = 5;
%M = 0.8;
uc=y(3);
i1=y(1);
uL1 = e - uc - R1*i1 - M_prev*dy2_prev; %takes M from previous step
M = fM(uL1);
D1 = L1/M-M/L2;
D2 = M/L1-L2/M;
%t_rect = mod(t,3);
%if t_rect < 1.5
% e = 120;
%else
% e = 0;
%end
%e=240*sin(t)
A = [-R1/(M*D1) R2/(L2*D1) -1/(M*D1)
-R1/(L1*D2) R2/(M*D2) -1/(L1*D2)
1/C 0 0];
b = [1/(M*D1)
1/(L1*D2)
0];
dy = A*y + b*e;
end
function M = fM(uL1)
global a;
[ud, md] = m_data();
if abs(uL1)>300
M=0.18;
else
M=polynewton(a, ud,abs(uL1));
end
end