    clear all
    clc
    close all
% 
% part1
% %create time vector interval
% t = (-5:0.001:5);
% 
% u1 = t>=0 %u(t)
% % 3 rows, 1 column and figure 1
% subplot(3,1,1)
% plot(t,u1)
% title('u(t)')
% xlabel('t')
% ylabel('u(t)')
% 
% u2 = t>=1 %u(t-1)
% % 3 rows, 1 comlumn and figure 2
% subplot(312)
% plot(t,u2)
% title('u(t-1)')
% xlabel('t')
% ylabel('u(t)')
% 
% u3 = t>=2 %u(t-2)
% % 3 rows, 1 comlumn and figure 3
% subplot(313)
% plot(t,u3)
% title('u(t-2)')
% xlabel('t')
% ylabel('u(t)')
% 
%Part 2
% %create time vector interval
% t = (-5:0.001:5);
% % u(t)−u(t−1)
% y1= (t>=0)-(t>=1);
% subplot(221)
% plot(t,y1)
% title('u(t)-u(t-1)')
% xlabel('t')
% ylabel('u(t)')
% 
% %u(t) +u(t−1) +u(t−2)
% y2 = (t>=0)+(t>=1)+(t>=2);
% subplot(222)
% plot(t,y2)
% title('u(t)+u(t-1)+u(t-2)')
% xlabel('t')
% ylabel('u(t)')
% 
% %u(t)u(t−2)
% y3 = (t>=0).*(t>=2);
% subplot(223)
% plot(t,y3)
% title('u(t)u(t-2)')
% xlabel('t')
% ylabel('u(t)')
% 
% %u(t) +u(t−1)−2u(t−2)
% y4 = (t>=0) + (t>=1) - 2*(t>=2);
% subplot(224)
% plot(t,y4)
% title('u(t)+u(t-1)-2u(t-2)')
% xlabel('t')
% ylabel('u(t)')
% 
% %part 3 
% t = (-2:0.01:2)
% v0 = (t>=0).*cos(2*pi*t);
% v1 = (t>=0).*cos(2*pi*t + pi/4);
% v2 = (t>=0).*cos(0.5*pi*t);
% v3 = v0.*v2
% v4 = v0.*(exp(-abs(t)))
% 
% subplot(2,2,1)
% plot(t,v0)
% hold on
% plot(t,v1)
% hold off
% grid on
% title('v0(t) vs v1(t)')
% xlabel('t')
% ylabel('u(t)')
% 
% subplot(222)
% plot(t,v0)
% hold on
% plot(t,v2)
% hold off
% grid on
% title('v0(t) vs v2(t)')
% xlabel('t')
% ylabel('u(t)')
% 
% subplot(223)
% plot(t,v0)
% hold on
% plot(t,v3)
% hold off
% grid on
% title('v0(t) vs v3(t)')
% xlabel('t')
% ylabel('u(t)')
% 
% subplot(224)
% plot (t,v0)
% hold on
% plot(t,v4)
% hold off
% grid on
% title('v0(t) vs v4(t)')
% xlabel('t')
% ylabel('u(t)')
% 
% part4
% s = tf('s')
% H1 = 1/(s+1);
%  H2 = (2*s^2+3*s+2)/(s^3+2*s+1);
% H = H1*H2
% 
% F1 = [1 5 6]
% poly2str(F1,'s')
% 
% s = 1
% polyval(F1,s)
% 
% roots(F1)
% 
% x= [ -1+2j -1-2j -1 0 0 ]
% F2 = poly(x)
% poly2str(F2,'s')
% % 
% H = conv(F1,F2)
%  poly2str(H,'s')
% 
% n = [1,1]
% d = [1,4,0]
% H3 = tf(n,d)
% % 
% [r,p,k] = residue(n,d)
% H3 = (r(1)/(s-p(1))) + (r(2)/(s-p(2)))
% syms s;
% ilaplace((s+1)/(s^2+4*s))
% 
% s = tf('s')
% num = [1,0,4]
% den = [1,4,0,0]
% H4 = tf(num,den)
% [r,p,k] = residue(num,den)
% syms s;
% ilaplace((s^2+4)/(s^3+4*s^2))

% n = [4,1,0,4]
% d=conv([1 0 0 0],conv ([1 3],conv ([1 0 2],[1 3]))); %Use "conv" to multiply polynomial
% H5 = tf(n,d)    
% [r,p,k] = residue(n,d)
% syms s;
% ilaplace((4*s^3+s^2+4)/(s^3*(s+3)^2*(s^2+2)))
% syms s            % create variable
% x=(4*s^3+s^2+4)/(s^3*(s+3)^2*(s^2+2)); % input equation
% pretty(x)         % to see equation in better way
% F=ilaplace(x);           % find inverse laplace
% pretty(F)  

% 
% s = tf('s');
% R1 = 50;
% L = 2*10^-3;
% R2 = 100;
% Req = (R1*R2)/(R1 + R2);
% H1 = s/(s+(R1/L));
% 
% %H2 = ((Req/R1)*s)/(s+Req/L);
% 
% 
% %hold on
% bode(H1);
% hold on
% %bode(H2)
% plot(25000,45,'*')
% 
% 
% legend('H1')
% title('Problem 5d')
w = 0:100000000;
R1 = 50;
L = 2*10^-3;
R2 = 100;
H1 = (1i*w)./((1i*w)+(R1/L));
H2 = (1i*w*2)./(
magH1 = abs(H1);
magH2 = abs(H2);
subplot(211)
semilogx(w,magH1);
hold on
xline(25000,'r')
hold on
yline(1, 'g')
hold on 
semilogx(w,magH2);
hold on 
xline((50000/3), 'r')
hold on 
yline((2/3), 'g')
legend( 'magH1', 'cut off freq', 'Hmax')
ylabel('magnitude')
xlabel('frequency')

Angle = angle(H1)*(180/pi);
subplot(212)
semilogx(w,Angle)
hold on
xline(25000,'r')
hold on
yline(45, 'g')
legend('angleH1','cut off freq', 'angleWc')
title('problem 5d')
ylabel('angle(jW)')
xlabel('frequency')


