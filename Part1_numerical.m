function Rlc
%h = 0.0001;
h=0.1
T = 0:h:30;
y = zeros(3, length(T));
e = zeros(1, length(T));
P= zeros(1, length(T));
dy = [0
0
0];
R1 = 0.1;
R2 = 10;
for i = 1:length(T) - 1
% E = 120*sin (2*pi*50*T(i)); %checking for different sin
% E = 210*sin (2*pi*5*T(i));
E=240*sin(T(i));
% E=1;
P(i+1)=R1*y(1,i)^2+R2*y(2,i)^2;
e(i) = E;
k1 = f(T(i), y(:,i),E);
k2 = f(T(i)+h/2, y(:,i) + h/2*k1,E);
k3 = f(T(i)+h/2, y(:,i) + h/2*k2,E);
k4 = f(T(i)+h, y(:,i) + h*k3,E);
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
grid;
end
function dy = f(t,y,e)
R1 = 0.1;
R2 = 10;
C = 0.5;
L1 = 3;
L2 = 5;
M = 0.8;
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