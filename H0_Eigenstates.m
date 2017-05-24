lamb=2.5;
E0 = 0.14285714285714286*lamb^(1/2);
E1 = 0.121207085423;
vec(1)=E0;
m1= [E0 E1
E1 0 ];
vec(2:3)=eig(m1);
vec(4:5)=eig(m1);
m2 = [E0 E1 E1 0
E1 0 0 E1
E1 0 0 E1
0 E1 E1 E0];
vec(6:9)=eig(m2);
vec(10:11)=[0 0];
vec(12:13)=eig(m1);
vec(14:15)=eig(m1);
vec(16)= E0;
Ens = vec - min(vec);