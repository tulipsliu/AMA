example7.mod_SparseAimDataR <- function(){

# example7.mod_aim_data()
#     This function will return various information about the AIM model,
#     but will not compute the G and H matrices.

  eqname = mat.or.vec(4, 1);
  param = mat.or.vec(5, 1);
  endog = mat.or.vec(4, 1);
  delay = mat.or.vec(4, 1);
  vtype = mat.or.vec(4, 1);
  eqtype = mat.or.vec(4, 1);

  modname = 'example7.mod';
  neq = 4;
  np = 5;
  nlag = 1;
  nlead = 1;

  eqname[1] = 'is';
  eqname[2] = 'phillips';
  eqname[3] = 'dr';
  eqname[4] = 'policy';

  eqtype[1] = 1;     eqtype[2] = 1;     eqtype[3] = 1;   
  eqtype[4] = 1;     eqtype_ = eqtype;

  param[1] = 'sigma';
  param[2] = 'delta';
  param[3] = 'lambda';
  param[4] = 'rho';
  param[5] = 'gampi';

  endog[1] = 'y';
  endog[2] = 'inf';
  endog[3] = 'dr';
  endog[4] = 'rs';

  delay[1] = 0;     delay[2] = 0;     delay[3] = 0;   
  delay[4] = 0;     delay_ = delay;

  vtype[1] = 1;     vtype[2] = 1;     vtype[3] = 1;   
  vtype[4] = 1;     vtype_ = vtype;

output <- list(eqname, param, endog, delay, vtype, eqtype)
output
}


