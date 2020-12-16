# example7.modAimMatricesR()
#     This script will compute the G and H matrices.

  h = mat.or.vec(4, 12);

  h[17] = h[17] + 1;
  h[29] = h[29] - (-1.0*((1.0*(sigma^-1.0))*1));
  h[33] = h[33] - 1;
  h[37] = h[37] - (-1.0*((1.0*(sigma^-1.0))*(-1.0*1)));
  h[22] = h[22] + 1;
  h[18] = h[18] - (lambda*1);
  h[38] = h[38] - (delta*1);
  h[27] = h[27] + 1;
  h[31] = h[31] - 1;
  h[15] = h[15] - (-1.0*1);
  h[32] = h[32] + 1;
  h[16] = h[16] - (rho*1);
  h[24] = h[24] - (gampi*1);

  cofh = h;
