ex_netcon = [1,2,3,4;
             4,3,2,1;
             5,6,7,8;
             8,7,6,5];
ex_gsyn = [0.2,0.2,0.2,0.2;
           0.2,0.2,0.2,0.2;
           0.2,0.2,0.2,0.2;
           0.2,0.2,0.2,0.2];

ex_flatg = [0.2,0.2,0.2,0.2];
s = 5;
psc = ex_gsyn .* (s * ex_netcon)
psc3 = s * (ex_netcon .* ex_gsyn)
rc = s * (ex_netcon .* ex_flatg)

