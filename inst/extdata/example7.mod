MODEL > EXAMPLE7

ENDOG> 

          y   _NOTD
        inf   _NOTD
         dr   _NOTD
         rs   _NOTD
  

EQUATION >  is
EQTYPE >    IMPOSED
EQ >        y = LEAD(y,1) - (1/sigma)*(rs-LEAD(inf,1)) 


EQUATION >  phillips
EQTYPE >    IMPOSED
EQ >        inf = delta*LEAD(inf,1) + lambda*y

EQUATION >  dr
EQTYPE >    IMPOSED
EQ >        dr = rs - LAG(rs,1)

EQUATION >  policy
EQTYPE >    IMPOSED
EQ >        rs  =  rho*LAG(rs,1) + gampi*inf

END
