>> E=1.4285714285714286

E =

    1.4286

>> v=[2 1 1 0]

v =

     2     1     1     0

>> M = E*diag(v)

M =

    2.8571         0         0         0
         0    1.4286         0         0
         0         0    1.4286         0
         0         0         0         0

>> E_2 = 0.19164522924380614

E_2 =

    0.1916

>>  v=[1 0]

v =

     1     0

>> M = E*diag(v)

M =

    1.4286         0
         0         0

>> E= 0.14285714285714288

E =

    0.1429

>> M = E*diag(v)

M =

    0.1429         0
         0         0

>> M=[M 0

0 0]
Error using horzcat
Dimensions of matrices being concatenated are not consistent.
 
>> M2=[M 0
0 M]
Error using horzcat
Dimensions of matrices being concatenated are not consistent.
 
>> zeros(2)

ans =

     0     0
     0     0

>> M2=[M zeros(2)
zeros(2) M]

M2 =

    0.1429         0         0         0
         0         0         0         0
         0         0    0.1429         0
         0         0         0         0

>> M3 = [ 0 1 0 0
0 0 0 0
1 0 0 1
0 1 0 0]

M3 =

     0     1     0     0
     0     0     0     0
     1     0     0     1
     0     1     0     0

>> M4 = M3 +M3.'

M4 =

     0     1     1     0
     1     0     0     1
     1     0     0     1
     0     1     1     0

>> FinM = M2+E_2*M4

FinM =

    0.1429    0.1916    0.1916         0
    0.1916         0         0    0.1916
    0.1916         0    0.1429    0.1916
         0    0.1916    0.1916         0

>> eig(FinM)

ans =

   -0.3247
    0.0586
    0.0843
    0.4676

>> E =1.4285714285714288

E =

    1.4286

>> E_2 = 1.9164522924380614

E_2 =

    1.9165

>> M = E*diag(v)

M =

    1.4286         0
         0         0

>> M2=[M zeros(2)
zeros(2) M]

M2 =

    1.4286         0         0         0
         0         0         0         0
         0         0    1.4286         0
         0         0         0         0

>> FinM = M2+E_2*M4

FinM =

    1.4286    1.9165    1.9165         0
    1.9165         0         0    1.9165
    1.9165         0    1.4286    1.9165
         0    1.9165    1.9165         0

>> eig(FinM)

ans =

   -3.2474
    0.5855
    0.8431
    4.6760

>> eig(FinM) = FinM/10
Subscript indices must either be real positive integers or logicals.
 
>> ens= eig(FinM)/10

ens =

   -0.3247
    0.0586
    0.0843
    0.4676

>> ens(1)

ans =

   -0.3247

>> ens= ens-ens(1)

ens =

         0
    0.3833
    0.4090
    0.7923

>> lambda = 2.5

lambda =

    2.5000

>> FinM = 2.5^(1/2)*M2+E_2*M4

FinM =

    2.2588    1.9165    1.9165         0
    1.9165         0         0    1.9165
    1.9165         0    2.2588    1.9165
         0    1.9165    1.9165         0

>> FinM = 2.5^(1/2)*M2 +E_2*M4

FinM =

    2.2588    1.9165    1.9165         0
    1.9165         0         0    1.9165
    1.9165         0    2.2588    1.9165
         0    1.9165    1.9165         0

>> ens= eig(FinM)/10

ens =

   -0.3012
    0.0821
    0.1437
    0.5270

>> ens= ens-ens(1)

ens =

         0
    0.3833
    0.4449
    0.8282

>> ens*0.5714285714285715

ans =

         0
    0.2190
    0.2542
    0.4732

>> 
