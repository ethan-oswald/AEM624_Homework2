clear;
clc;

M1 = 3;
gamma = 1.4;
P1 = 101325;
T1 = 298.15;
x(1) = .01370207;
y(1) = 0;
dy(1) = 0;
q1 = gamma/2 * P1 * M1^2;

for i = 2:241
x(i) = i * 0.01370207;
y(i) = -0.008333 + 0.609425*x(i) - 0.092593*x(i)^2;
dy(i) = 0.609425 - 0.185186* x(i);
theta(i) = atand(dy(i));
end

theta_n = theta(2);

testtheta = -999;
testB = 0;

while (abs(testtheta-theta_n)) > 0.01
testtheta = atand(2*cotd(testB) * ((M1^2*sind(testB)^2-1)/(M1^2*(gamma+cosd(2*testB))+2)));
testB = testB + .0001;
end
Beta = testB;
Machn1 = M1 * sind(Beta);
Mn2 = sqrt((Machn1^2 + (2/(gamma-1)))/((2*gamma/(gamma-1)*Machn1^2)-1));
Mn = Mn2 / (sind(Beta - theta_n));
nu1 = sqrt((gamma+1)/(gamma-1))*atand(sqrt((gamma-1)/(gamma+1)*(Mn^2-1)))-atand(sqrt(Mn^2-1));
pn_ratio = 1 + 2 * gamma/(gamma+1) * (Machn1^2-1);


for i = 2:241

%Newtonian Method
c_p_newt(i) = 2 * sind(theta(i))^2;
p2_newt(i) = q1 * c_p_newt(i) + P1;
p2_p1_newtonian(i) = p2_newt(i) / P1;

%%Modified Newtonian Method
c_p_max = ((gamma+1)^2/(4*gamma))^(gamma/(gamma-1))*4/(gamma+1);
c_p_mod(i) = c_p_max * sind(theta(i))^2;
p2_mod(i) = q1 * c_p_mod(i) + P1;
p2_p1_modified(i) = p2_mod(i) / P1;

%%Tangent Wedge Method
testtheta = -999;
testB = 0;

while (abs(testtheta-theta(i))) > 0.01
testtheta = atand(2*cotd(testB) * ((M1^2*sind(testB)^2-1)/(M1^2*(gamma+cosd(2*testB))+2)));
testB = testB + .0001;
end
B(i) = testB;
Mn1(i) = M1 * sind(B(i));
p2_p1_TW(i) = 1 + 2*gamma/(gamma+1) * (Mn1(i)^2-1);

%%Shock-Expansion Method

testM = 0.1;
testnu = -999;
nu2(i) = -theta(i) + nu1 + theta_n;
while (abs(testnu-nu2(i))) > 0.01
testM = testM + .0001;
testnu = sqrt((gamma+1)/(gamma-1))*atand(sqrt(((gamma-1)/(gamma+1))*(testM^2-1)))-atand(sqrt(testM^2-1));
end
Mi(i) = testM;

pi_pn(i) = ((1+(gamma-1)/2 * Mn^2)/(1+(gamma-1)/2 * Mi(i)^2))^(gamma/(gamma-1));

p2_p1_SE(i) = pi_pn(i) * pn_ratio;
end

plot(x, p2_p1_newtonian)
hold on 
plot(x, p2_p1_modified)
plot(x, p2_p1_TW)
plot(x,p2_p1_SE)
axis([.028 3.5 0 9]);
legend("Newtonian Method", "Modified Newtonian Method", "Tangent Wedge Method", "Shock-Expansion Method")
title("Pressure Ratio v. Distance along the Body")
ylabel("P2/P1")
xlabel("Distance along the body")
hold off



