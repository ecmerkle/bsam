/* This file is based on LERSIL.stan by Ben Goodrich.
   https://github.com/bgoodri/LERSIL                  */
functions { // you can use these in R following `rstan::expose_stan_functions("foo.stan")`
  /*
    Fills in the elements of a coefficient matrix containing some mix of 
    totally free, free subject to a sign constraint, and fixed elements
    
    @param free_elements vector of unconstrained elements
    @param skeleton matrix of the same dimensions as the output whose elements are
      positive_infinity(): if output element is totally free
      other: if output element is fixed to that number
    @return matrix of coefficients
  */
  matrix fill_matrix(vector free_elements, matrix skeleton, int[,] eq_skeleton, int pos_start) {
    int R = rows(skeleton);
    int C = cols(skeleton);
    matrix[R, C] out;

    int pos = pos_start;
    int eqelem = 0;
    for (c in 1:C) for (r in 1:R) {
      real rc = skeleton[r, c];
      if (is_inf(rc)) { // free
	real eq = eq_skeleton[pos, 1];
	if (eq == 0) {
	  out[r,c] = free_elements[pos];
	} else {
	  eqelem = eq_skeleton[pos, 2];
	  out[r,c] = free_elements[eqelem];
	}
	pos += 1;
      } else out[r,c] = skeleton[r, c]; // fixed, so do not bump pos
    }
    return out;
  }
  
  /*
   * This is a bug-free version of csr_to_dense_matrix and has the same arguments
   */
  matrix to_dense_matrix(int m, int n, vector w, int[] v, int[] u) {
    matrix[m, n] out = rep_matrix(0, m, n);
    int pos = 1;
    for (i in 1:m) {
      int start = u[i];
      int nnz = u[i + 1] - start;
      for (j in 1:nnz) {
        out[i, v[pos]] = w[pos];
        pos += 1;
      }
    }
    return out;
  }

  // sign function
  int sign(real x) {
    if (x > 0)
      return 1;
    else
      return -1;
  }

  // sign-constrain a vector of loadings
  vector sign_constrain_load(vector free_elements, int npar, int[,] sign_mat) {
    vector[npar] out;
    for (i in 1:npar) {
      if (sign_mat[i,1]) {
        int lookupval = sign_mat[i,2];
        if (free_elements[lookupval] < 0) {
	  out[i] = -free_elements[i];
	} else {
	  out[i] = free_elements[i];
	}
      } else {
        out[i] = free_elements[i];
      }
    }
    return out;
  }

  // sign-constrain a vector of regressions or covariances
  vector sign_constrain_reg(vector free_elements, int npar, int[,] sign_mat, vector load_par1, vector load_par2) {
    vector[npar] out;
    for (i in 1:npar) {
      if (sign_mat[i,1]) {
        int lookupval1 = sign_mat[i,2];
	int lookupval2 = sign_mat[i,3];
        if (sign(load_par1[lookupval1]) * sign(load_par2[lookupval2]) < 0) {
	  out[i] = -free_elements[i];
	} else {
	  out[i] = free_elements[i];
	}
      } else {
        out[i] = free_elements[i];
      }
    }
    return out;
  }

  // obtain covariance parameter vector for correlation/sd matrices
  vector cor2cov(matrix[] cormat, matrix[] sdmat, vector free_elements, int[,] wskel, int ngrp) {
    vector[num_elements(free_elements)] out;
    int R = rows(to_matrix(cormat[1]));
    int C = cols(to_matrix(cormat[1]));
    int pos = 1;
    for (g in 1:ngrp) {
      for (c in 1:(R-1)) for (r in (c+1):R) {
        if (cormat[g,r,c] != 0 && wskel[pos,1] == 0) {
	  out[pos] = sdmat[g,r,r] * sdmat[g,c,c] * cormat[g,r,c];
	  pos += 1;
	}
      }
    }
    return out;
  }
  
}
data {
  // see https://books.google.com/books?id=9AC-s50RjacC&lpg=PP1&dq=LISREL&pg=PA2#v=onepage&q=LISREL&f=false
  int<lower=0> p; // number of manifest response variables
  int<lower=0> q; // number of manifest predictors
  int<lower=0> m; // number of latent endogenous variables
  int<lower=0> n; // number of latent exogenous variables
  int<lower=1> Ng; // number of groups
  cov_matrix[p + q] S[Ng];     // sample covariance matrix among all manifest variables NB!! multiply by (N-1) to use wishart lpdf!!
  int<lower=0, upper=1> has_data; // are the raw (centered) data on y and x available?
  int<lower=0, upper=1> missing; // are there missing values?
  int<lower=0, upper=1> save_lvs; // should we save lvs?
  int<lower=1> Np; // number of group-by-missing patterns combos
  int<lower= 1> N[Ng]; // number of observations per group
  int<lower=1> Nobs[Np]; // number of observed variables in each missing pattern
  int<lower=0> Obsvar[Np, p + q];
  int<lower=1> Ntot; // number of observations across all groups
  int<lower=1> startrow[Np]; // starting row for each missing pattern
  int<lower=1,upper=Ntot> endrow[Np]; // ending row for each missing pattern
  int<lower=1,upper=Ng> grpnum[Np]; // group number for each row of data
  vector[p + q] YX[has_data ? Ntot : 0]; // if data, include them

  
  /* sparse matrix representations of skeletons of coefficient matrices, 
     which is not that interesting but necessary because you cannot pass
     missing values into the data block of a Stan program from R */
  int<lower=0> len_w1;        // number of free elements in Lambda_y
  vector[len_w1] w1;          // values of free elements in Lambda_y
  int<lower=1> v1[len_w1];    // index  of free elements in Lambda_y
  int<lower=1> u1[p + 1];     // index  of free elements in Lambda_y
  int<lower=0> w1skel[Ng * len_w1, 2];
  int<lower=0> lam_y_sign[Ng * len_w1, 2];
  int<lower=0> len_lam_y;     // number of free elements minus equality constraints
  real lambda_y_mn[len_lam_y];           // prior
  real<lower=0> lambda_y_sd[len_lam_y];

  // same things but for Lambda_x
  int<lower=0> len_w2;
  vector[len_w2] w2;
  int<lower=1> v2[len_w2];
  int<lower=1> u2[q + 1];
  int<lower=0> w2skel[Ng * len_w2, 2];
  int<lower=0> lam_x_sign[Ng * len_w2, 2];
  int<lower=0> len_lam_x;
  real lambda_x_mn[len_lam_x];
  real<lower=0> lambda_x_sd[len_lam_x];
  
  // same things but for Gamma
  int<lower=0> len_w3;
  vector[len_w3] w3;
  int<lower=1> v3[len_w3];
  int<lower=1> u3[m + 1];
  int<lower=0> w3skel[Ng * len_w3, 2];
  int<lower=0> gam_sign[Ng * len_w3, 3];
  int<lower=0> len_gam;
  real gamma_mn[len_gam];
  real<lower=0> gamma_sd[len_gam];
  
  // same things but for B
  int<lower=0> len_w4;
  vector[len_w4] w4;
  int<lower=1> v4[len_w4];
  int<lower=1> u4[m + 1];
  int<lower=0> w4skel[Ng * len_w4, 2];
  int<lower=0> b_sign[Ng * len_w4, 3];
  int<lower=0> len_b;
  real b_mn[len_b];
  real<lower=0> b_sd[len_b];
  
  // same things but for diag(Theta)
  int<lower=0> len_w5;
  vector[len_w5] w5;
  int<lower=1> v5[len_w5];
  int<lower=1> u5[p + 1];
  int<lower=0> w5skel[Ng * len_w5, 2];
  int<lower=0> len_thet_sd;
  real theta_sd_rate[len_thet_sd];

  // same things but for diag(Theta_x)
  int<lower=0> len_w6;
  vector[len_w6] w6;
  int<lower=1> v6[len_w6];
  int<lower=1> u6[q + 1];
  int<lower=0> w6skel[Ng * len_w6, 2];
  int<lower=0> len_thet_x_sd;
  real theta_x_sd_rate[len_thet_x_sd];
  
  // same things but for Theta_r
  int<lower=0> len_w7;
  vector[len_w7] w7;
  int<lower=1> v7[len_w7];
  int<lower=1> u7[p + 1];
  int<lower=0> w7skel[Ng * len_w7, 2];
  int<lower=0> len_thet_r;
  real<lower=0> theta_r_alpha[len_thet_r];
  real<lower=0> theta_r_beta[len_thet_r];
  
  // same things but for Theta_r_x
  int<lower=0> len_w8;
  vector[len_w8] w8;
  int<lower=1> v8[len_w8];
  int<lower=1> u8[q + 1];
  int<lower=0> w8skel[Ng * len_w8, 2];
  int<lower=0> len_thet_x_r;
  real<lower=0> theta_x_r_alpha[len_thet_x_r];
  real<lower=0> theta_x_r_beta[len_thet_x_r];
  
  // same things but for Psi
  int<lower=0> len_w9;
  vector[len_w9] w9;
  int<lower=1> v9[len_w9];
  int<lower=1> u9[m + 1];
  int<lower=0> w9skel[Ng * len_w9, 2];
  int<lower=0> len_psi_sd;
  real<lower=0> psi_sd_rate[len_psi_sd];

  // same things but for Psi_r
  int<lower=0> len_w10;
  vector[len_w10] w10;
  int<lower=1> v10[len_w10];
  int<lower=1> u10[m + 1];
  int<lower=0> w10skel[Ng * len_w10, 2];
  int<lower=0> psi_r_sign[Ng * len_w10, 3];
  int<lower=0> len_psi_r;
  real<lower=0> psi_r_alpha[len_psi_r];
  real<lower=0> psi_r_beta[len_psi_r];
  
  // same things but for Phi
  int<lower=0> len_w11;
  vector[len_w11] w11;
  int<lower=1> v11[len_w11];
  int<lower=1> u11[n + 1];
  int<lower=0> w11skel[Ng * len_w11, 2];
  int<lower=0> len_phi_sd;
  real<lower=0> phi_sd_rate[len_phi_sd];

  // same things but for Phi_r
  int<lower=0> len_w12;
  vector[len_w12] w12;
  int<lower=1> v12[len_w12];
  int<lower=1> u12[n + 1];
  int<lower=0> w12skel[Ng * len_w12, 2];
  int<lower=0> phi_r_sign[Ng * len_w12, 3];
  int<lower=0> len_phi_r;
  real<lower=0> phi_r_alpha[len_phi_r];
  real<lower=0> phi_r_beta[len_phi_r];
  
  // same things but for Nu
  int<lower=0> len_w13;
  vector[len_w13] w13;
  int<lower=1> v13[len_w13];
  int<lower=1> u13[p + q + 1];
  int<lower=0> w13skel[Ng * len_w13, 2];
  int<lower=0> len_nu;
  real nu_mn[len_nu];
  real<lower=0> nu_sd[len_nu];
  
  // same things but for Alpha
  int<lower=0> len_w14;
  vector[len_w14] w14;
  int<lower=1> v14[len_w14];
  int<lower=1> u14[m + n + 1];
  int<lower=0> w14skel[Ng * len_w14, 2];
  int<lower=0> len_alph;
  real alpha_mn[len_alph];
  real<lower=0> alpha_sd[len_alph];
}
transformed data { // (re)construct skeleton matrices in Stan (not that interesting)
  matrix[p, m] Lambda_y_skeleton = to_dense_matrix(p, m, w1, v1, u1);
  matrix[q, n] Lambda_x_skeleton = to_dense_matrix(q, n, w2, v2, u2);
  matrix[m, n] Gamma_skeleton = to_dense_matrix(m, n, w3, v3, u3);
  matrix[m, m] B_skeleton = to_dense_matrix(m, m, w4, v4, u4);
  matrix[p, p] Theta_skeleton = to_dense_matrix(p, p, w5, v5, u5);
  matrix[q, q] Theta_x_skeleton = to_dense_matrix(q, q, w6, v6, u6);
  matrix[p, p] Theta_r_skeleton = to_dense_matrix(p, p, w7, v7, u7);
  matrix[q, q] Theta_x_r_skeleton = to_dense_matrix(q, q, w8, v8, u8);
  matrix[m, m] Psi_skeleton = to_dense_matrix(m, m, w9, v9, u9);
  matrix[m, m] Psi_r_skeleton = to_dense_matrix(m, m, w10, v10, u10);
  matrix[n, n] Phi_skeleton = to_dense_matrix(n, n, w11, v11, u11);
  matrix[n, n] Phi_r_skeleton = to_dense_matrix(n, n, w12, v12, u12);
  matrix[(p + q), 1] Nu_skeleton = to_dense_matrix((p + q), 1, w13, v13, u13);
  matrix[(m + n), 1] Alpha_skeleton = to_dense_matrix((m + n), 1, w14, v14, u14);

  matrix[m, m] I = diag_matrix(rep_vector(1, m));

  int len_free1 = 0;
  int len_free2 = 0;
  int len_free3 = 0;
  int len_free4 = 0;
  int len_free5 = 0;
  int len_free6 = 0;
  int len_free7 = 0;
  int len_free8 = 0;
  int len_free9 = 0;
  int len_free10 = 0;
  int len_free11 = 0;
  int len_free12 = 0;
  int len_free13 = 0;
  int len_free14 = 0;

  int g_start1[Ng];
  int g_start2[Ng];
  int g_start3[Ng];
  int g_start4[Ng];
  int g_start5[Ng];
  int g_start6[Ng];
  int g_start7[Ng];
  int g_start8[Ng];
  int g_start9[Ng];
  int g_start10[Ng];
  int g_start11[Ng];
  int g_start12[Ng];
  int g_start13[Ng];
  int g_start14[Ng];

  int pos;
  
  for (g in 1:Ng) {
    // count free elements in Lambda_y_skeleton
    pos = len_free1 + 1;
    g_start1[g] = pos;
    for (i in 1:p) {
      for (j in 1:m) {
        if (is_inf(Lambda_y_skeleton[i,j])) {
	  if (w1skel[pos,1] == 0) len_free1 += 1;
	  pos += 1;
        }
      }
    }

    // same thing but for Lambda_x_skeleton
    pos = len_free2 + 1;
    g_start2[g] = pos;
    for (i in 1:q) {
      for (j in 1:n) {
	if (is_inf(Lambda_x_skeleton[i,j])) {
	  if (w2skel[pos,2] == 0) len_free2 += 1;
	  pos += 1;
	}
      }
    }
  
    // same thing but for Gamma_skeleton
    pos = len_free3 + 1;
    g_start3[g] = pos;
    for (i in 1:m) {
      for (j in 1:n) {
	if (is_inf(Gamma_skeleton[i,j])) {
	  if (w3skel[pos,2] == 0) len_free3 += 1;
	  pos += 1;
	}
      }
    }

    // same thing but for B_skeleton
    pos = len_free4 + 1;
    g_start4[g] = pos;
    for (i in 1:m) {
      for (j in 1:m) {
	if (is_inf(B_skeleton[i,j])) {
	  if (w4skel[pos,2] == 0) len_free4 += 1;
	  pos += 1;
	}
      }
    }
    
    // same thing but for Theta_skeleton
    pos = len_free5 + 1;
    g_start5[g] = pos;
    for (i in 1:p) {
      if (is_inf(Theta_skeleton[i,i])) {
	if (w5skel[pos,2] == 0) len_free5 += 1;
	pos += 1;
      }
    }

    // same thing but for Theta_x_skeleton
    pos = len_free6 + 1;
    g_start6[g] = pos;
    for (i in 1:q) {
      if (is_inf(Theta_x_skeleton[i,i])) {
	if (w6skel[pos,2] == 0) len_free6 += 1;
	pos += 1;
      }
    }

    // same thing but for Theta_r_skeleton
    pos = len_free7 + 1;
    g_start7[g] = pos;
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
	if (is_inf(Theta_r_skeleton[j,i])) {
	  if (w7skel[pos,2] == 0) len_free7 += 1;
	  pos += 1;
	}
      }
    }

    // same thing but for Theta_x_r_skeleton
    pos = len_free8 + 1;
    g_start8[g] = pos;
    for (i in 1:(q-1)) {
      for (j in (i+1):q) {
	if (is_inf(Theta_x_r_skeleton[j,i])) {
	  if (w8skel[pos,2] == 0) len_free8 += 1;
	  pos += 1;
	}
      }
    }

    // same thing but for Psi_skeleton
    pos = len_free9 + 1;
    g_start9[g] = pos;
    for (i in 1:m) {
      if (is_inf(Psi_skeleton[i,i])) {
	if (w9skel[pos,2] == 0) len_free9 += 1;
	pos += 1;
      }
    }

    // same thing but for Psi_r_skeleton
    pos = len_free10 + 1;
    g_start10[g] = pos;
    for (i in 1:(m-1)) {
      for (j in (i+1):m) {
	if (is_inf(Psi_r_skeleton[j,i])) {
	  if (w10skel[pos,2] == 0) len_free10 += 1;
	  pos += 1;
	}
      }
    }

    // same thing but for Phi_skeleton
    pos = len_free11 + 1;
    g_start11[g] = pos;
    for (i in 1:n) {
      if (is_inf(Phi_skeleton[i,i])) {
	if (w11skel[pos,2] == 0) len_free11 += 1;
	pos += 1;
      }
    }

    // same thing but for Phi_r_skeleton
    pos = len_free12 + 1;
    g_start12[g] = pos;
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
	if (is_inf(Phi_r_skeleton[j,i])) {
	  if (w12skel[pos,2] == 0) len_free12 += 1;
	  pos += 1;
	}
      }
    }

    // same thing but for Nu_skeleton
    pos = len_free13 + 1;
    g_start13[g] = pos;
    for (i in 1:(p+q)) {
      if (is_inf(Nu_skeleton[i,1])) {
	if (w13skel[pos,2] == 0) len_free13 += 1;
	pos += 1;
      }
    }

    // same thing but for Alpha_skeleton
    pos = len_free14 + 1;
    g_start14[g] = pos;
    for (i in 1:(m+n)) {
      if (is_inf(Alpha_skeleton[i,1])) {
	if (w14skel[pos,2] == 0) len_free14 += 1;
	pos += 1;
      }
    }
  }
}
parameters {
  // free elements (possibly with inequality constraints) for coefficient matrices
  vector[len_free1] Lambda_y_free;
  vector[len_free2] Lambda_x_free;
  vector[len_free3] Gamma_free;
  vector[len_free4] B_free;
  vector<lower=0>[len_free5] Theta_sd_free;
  vector<lower=0>[len_free6] Theta_x_sd_free;
  vector<lower=0,upper=1>[len_free7] Theta_r_free; // to use beta prior
  vector<lower=0,upper=1>[len_free8] Theta_x_r_free;
  vector<lower=0>[len_free9] Psi_sd_free;
  vector<lower=0,upper=1>[len_free10] Psi_r_free;
  vector<lower=0>[len_free11] Phi_sd_free;
  vector<lower=0,upper=1>[len_free12] Phi_r_free;
  vector[len_free13] Nu_free;
  vector[len_free14] Alpha_free;
}
transformed parameters {
  matrix[p, m] Lambda_y[Ng];
  matrix[q, n] Lambda_x[Ng];
  matrix[m, n] Gamma[Ng];
  matrix[m, m] B[Ng];
  matrix[p, p] Theta_sd[Ng];
  matrix[q, q] Theta_x_sd[Ng];
  matrix[p, p] T_r_lower[Ng];
  matrix[p, p] Theta_r[Ng];
  matrix[q, q] T_x_r_lower[Ng];
  matrix[q, q] Theta_x_r[Ng];
  matrix[p + q, 1] Nu[Ng];
  matrix[m + n, 1] Alpha[Ng];

  matrix[m, m] Psi[Ng];
  matrix[n, n] PHI[Ng];
  
  matrix[m, m] Psi_sd[Ng];
  matrix[m, m] Psi_r_lower[Ng];
  matrix[m, m] Psi_r[Ng];
  matrix[n, n] Phi_sd[Ng];
  matrix[n, n] Phi_r_lower[Ng];
  matrix[n, n] Phi_r[Ng];

  // Now fill them in
  for (g in 1:Ng) {
    Lambda_y[g] = fill_matrix(Lambda_y_free, Lambda_y_skeleton, w1skel, g_start1[g]);
    Lambda_x[g] = fill_matrix(Lambda_x_free, Lambda_x_skeleton, w2skel, g_start2[g]);
    Gamma[g] = fill_matrix(Gamma_free, Gamma_skeleton, w3skel, g_start3[g]);
    B[g] = fill_matrix(B_free, B_skeleton, w4skel, g_start4[g]);
    Theta_sd[g] = fill_matrix(Theta_sd_free, Theta_skeleton, w5skel, g_start5[g]);
    Theta_x_sd[g] = fill_matrix(Theta_x_sd_free, Theta_x_skeleton, w6skel, g_start6[g]);
    T_r_lower[g] = fill_matrix(2*Theta_r_free - 1, Theta_r_skeleton, w7skel, g_start7[g]);
    Theta_r[g] = T_r_lower[g] + transpose(T_r_lower[g]) - diag_matrix(rep_vector(1, p));
    T_x_r_lower[g] = fill_matrix(2*Theta_x_r_free - 1, Theta_x_r_skeleton, w8skel, g_start8[g]);
    Theta_x_r[g] = T_x_r_lower[g] + transpose(T_x_r_lower[g]) - diag_matrix(rep_vector(1, q));
    Nu[g] = fill_matrix(Nu_free, Nu_skeleton, w13skel, g_start13[g]);
    Alpha[g] = fill_matrix(Alpha_free, Alpha_skeleton, w14skel, g_start14[g]);

    Psi[g] = diag_matrix(rep_vector(0, m));
    PHI[g] = diag_matrix(rep_vector(0, n));
  
    if (m > 0) {
      Psi_sd[g] = fill_matrix(Psi_sd_free, Psi_skeleton, w9skel, g_start9[g]);
      Psi_r_lower[g] = fill_matrix(2*Psi_r_free - 1, Psi_r_skeleton, w10skel, g_start10[g]);
      Psi_r[g] = Psi_r_lower[g] + transpose(Psi_r_lower[g]) - diag_matrix(rep_vector(1, m));
      Psi[g] = quad_form_sym(Psi_r[g], Psi_sd[g]);
    }

    if (n > 0) {
      Phi_sd[g] = fill_matrix(Phi_sd_free, Phi_skeleton, w11skel, g_start11[g]);
      Phi_r_lower[g] = fill_matrix(2*Phi_r_free - 1, Phi_r_skeleton, w12skel, g_start12[g]);
      Phi_r[g] = Phi_r_lower[g] + transpose(Phi_r_lower[g]) - diag_matrix(rep_vector(1, n));
      PHI[g] = quad_form_sym(Phi_r[g], Phi_sd[g]);
    }
  }
}
model { // N.B.: things declared in the model block do not get saved in the output, which is okay here
  matrix[p, m] Lambda_y_A[Ng];     // = Lambda_y * (I - B)^{-1}
  matrix[m, m] GPG[Ng];
  matrix[n, q] Lambda_xt[Ng];                         // copies so do it just once

  // see https://books.google.com/books?id=9AC-s50RjacC&lpg=PP1&dq=LISREL&pg=PA3#v=onepage&q=LISREL&f=false
  vector[p + q] Mu[Ng];
  matrix[p + q, p + q] Sigma[Ng];                                           // model covariance matrix

  matrix[p, q] top_right[Ng];        // top right block of Sigma
  
  for (g in 1:Ng) {
    Lambda_y_A[g] = mdivide_right(Lambda_y[g], I - B[g]);     // = Lambda_y * (I - B)^{-1}
    Lambda_xt[g] = transpose(Lambda_x[g]);                         // copies so do it just once

    Mu[g] = to_vector(Nu[g]);

    if (q > 0) {
      top_right[g] = Lambda_y_A[g] * Gamma[g] * PHI[g] * Lambda_xt[g];        // top right block of Sigma    
      Sigma[g, 1:p, (p + 1):(p + q)] = top_right[g];
      Sigma[g, (p + 1):(p + q), 1:p] = transpose(top_right[g]);
      Sigma[g, (p + 1):(p + q), (p + 1):(p + q)] = quad_form(PHI[g], Lambda_xt[g]);
      Sigma[g, (p + 1):(p + q), (p + 1):(p + q)] += quad_form_sym(Theta_x_r[g], Theta_x_sd[g]);

      if (n > 0) {
	Mu[g, (p + 1):(p + q)] += to_vector(Lambda_x[g] * Alpha[g, (m + 1):(m + n), 1]);
      }
    }

    GPG[g] = diag_matrix(rep_vector(0, m));
    if (p > 0) {
      if (q > 0) {
	GPG[g] = quad_form(PHI[g], transpose(Gamma[g]));
      }
      Sigma[g, 1:p, 1:p] = quad_form(GPG[g] + Psi[g], transpose(Lambda_y_A[g]));
      Sigma[g, 1:p, 1:p] += quad_form_sym(Theta_r[g], Theta_sd[g]);

      if (m > 0) {
	Mu[g, 1:p] += to_vector(Lambda_y_A[g] * Alpha[g, 1:m, 1]);
      }
      if (n > 0) {
	Mu[g, 1:p] += to_vector(Lambda_y_A[g] * Gamma[g] * Alpha[g, (m + 1):(m + n), 1]);
      }
    }
  }
    
  /* log-likelihood */
  if (has_data) {
    int obsidx[p + q];
    int r1;
    int r2;
    int grpidx;
    for (mm in 1:Np) {
      obsidx = Obsvar[mm,];
      r1 = startrow[mm];
      r2 = endrow[mm];
      grpidx = grpnum[mm];
      target += multi_normal_lpdf(YX[r1:r2,1:Nobs[mm]] | Mu[grpidx, obsidx[1:Nobs[mm]]], Sigma[grpidx, obsidx[1:Nobs[mm]], obsidx[1:Nobs[mm]]]);
    }
  } else {
    for (g in 1:Ng) {
      target += wishart_lpdf(S[g] | N[g] - 1, Sigma[g]);
    }
  }
  
  /* prior densities in log-units */
  target += normal_lpdf(Lambda_y_free | lambda_y_mn, lambda_y_sd);
  target += normal_lpdf(Lambda_x_free | lambda_x_mn, lambda_x_sd);
  target += normal_lpdf(Gamma_free    | gamma_mn, gamma_sd);
  target += normal_lpdf(B_free        | b_mn, b_sd);

  target += normal_lpdf(Nu_free       | nu_mn, nu_sd);
  target += normal_lpdf(Alpha_free    | alpha_mn, alpha_sd);
  
  target += exponential_lpdf(Theta_sd_free | theta_sd_rate);
  target += exponential_lpdf(Theta_x_sd_free | theta_x_sd_rate);
  target += exponential_lpdf(Psi_sd_free | psi_sd_rate);
  target += exponential_lpdf(Phi_sd_free | phi_sd_rate);

  target += beta_lpdf(Theta_r_free | theta_r_alpha, theta_r_beta);
  target += beta_lpdf(Theta_x_r_free | theta_x_r_alpha, theta_x_r_beta);
  target += beta_lpdf(Psi_r_free | psi_r_alpha, psi_r_beta);
  target += beta_lpdf(Phi_r_free | phi_r_alpha, phi_r_beta);
}
generated quantities { // these matrices are saved in the output but do not figure into the likelihood
  // see https://books.google.com/books?id=9AC-s50RjacC&lpg=PP1&dq=LISREL&pg=PA34#v=onepage&q=LISREL&f=false

  matrix[Ntot, save_lvs ? m + n : 0] eta;
  // matrix[Ntot, has_data ? m : 0] eta;
  // matrix[Ntot, has_data ? n : 0] xi;

  // sign constraints and correlations
  vector[len_free1] ly_sign;
  matrix[p, m] L_Y[Ng];
  vector[len_free2] lx_sign;
  matrix[q, n] L_X[Ng];
  vector[len_free3] g_sign;
  matrix[m, n] Gam[Ng];
  vector[len_free4] bet_sign;
  matrix[m, m] Bet[Ng];
  matrix[p, p] Theta[Ng];
  matrix[q, q] Theta_x[Ng];
  matrix[m, m] PSmat[Ng];
  matrix[m, m] PS[Ng];
  matrix[n, n] PHmat[Ng];
  matrix[n, n] PH[Ng];
  vector[len_free7] Theta_cov;
  vector[len_free5] Theta_var;
  vector[len_free8] Theta_x_cov;
  vector[len_free6] Theta_x_var;
  vector[len_free10] P_r;
  vector[len_free10] Psi_cov;
  vector[len_free9] Psi_var;
  vector[len_free12] Ph_r;
  vector[len_free12] Ph_cov;
  vector[len_free11] Ph_var;

  // first deal with sign constraints:
  ly_sign = sign_constrain_load(Lambda_y_free, len_free1, lam_y_sign);
  lx_sign = sign_constrain_load(Lambda_x_free, len_free2, lam_x_sign);
  g_sign = sign_constrain_reg(Gamma_free, len_free3, gam_sign, Lambda_x_free, Lambda_y_free);
  bet_sign = sign_constrain_reg(B_free, len_free4, b_sign, Lambda_y_free, Lambda_y_free);
  P_r = sign_constrain_reg(2 * Psi_r_free - 1, len_free10, psi_r_sign, Lambda_y_free, Lambda_y_free);
  Ph_r = sign_constrain_reg(2 * Phi_r_free - 1, len_free12, phi_r_sign, Lambda_x_free, Lambda_x_free);
  
  for (g in 1:Ng) {
    L_Y[g] = fill_matrix(ly_sign, Lambda_y_skeleton, w1skel, g_start1[g]);

    L_X[g] = fill_matrix(lx_sign, Lambda_x_skeleton, w2skel, g_start2[g]);

    Gam[g] = fill_matrix(g_sign, Gamma_skeleton, w3skel, g_start3[g]);

    Bet[g] = fill_matrix(bet_sign, B_skeleton, w4skel, g_start4[g]);

    Theta[g] = quad_form_sym(Theta_r[g], Theta_sd[g]);

    if (q > 0) {
      Theta_x[g] = quad_form_sym(Theta_x_r[g], Theta_x_sd[g]);
    }

    if (m > 0) {
      PSmat[g] = fill_matrix(P_r, Psi_r_skeleton, w10skel, g_start10[g]);
      PS[g] = quad_form_sym(PSmat[g] + transpose(PSmat[g]) - diag_matrix(rep_vector(1, m)), Psi_sd[g]);
    }

    if (n > 0) {
      PHmat[g] = fill_matrix(Ph_r, Phi_r_skeleton, w12skel, g_start12[g]);
      PH[g] = quad_form_sym(PHmat[g] + transpose(PHmat[g]) - diag_matrix(rep_vector(1, n)), Phi_sd[g]);
    }
    
  }

  // off-diagonal covariance parameter vectors, from cor/sd matrices:
  Theta_cov = cor2cov(Theta_r, Theta_sd, Theta_r_free, w7skel, Ng);
  Theta_var = Theta_sd_free .* Theta_sd_free;
  Theta_x_cov = cor2cov(Theta_x_r, Theta_x_sd, Theta_x_r_free, w8skel, Ng);
  Theta_x_var = Theta_x_sd_free .* Theta_x_sd_free;
  if (m > 0) {
    /* iden is created so that we can re-use cor2cov, even though
       we don't need to multiply to get covariances */
    matrix[m, m] iden[Ng];
    for (g in 1:Ng) {
      iden[g] = diag_matrix(rep_vector(1, m));
    }
    Psi_cov = cor2cov(PS, iden, P_r, w10skel, Ng);
  }
  if (n > 0) {
    matrix[n, n] iden[Ng];
    for (g in 1:Ng) {
      iden[g] = diag_matrix(rep_vector(1, n));
    }
    Ph_cov = cor2cov(PH, iden, Ph_r, w12skel, Ng);
  }
  Psi_var = Psi_sd_free .* Psi_sd_free;
  Ph_var = Phi_sd_free .* Phi_sd_free;

  // now use matrices with sign fixes to deal with lvs
  if (save_lvs) { // all matrices defined in this local block are not saved in the output
    matrix[m, m] A;
    matrix[m, n] total_xi_eta;
    matrix[m, n] indirect_xi_eta;
    matrix[m, m] total_eta_eta;
    matrix[m, m] indirect_eta_eta;
    matrix[p, m] total_eta_y;
    matrix[p, m] indirect_eta_y;
    matrix[p, n] total_xi_y; // = indirect_xi_y since there is no direct effect

    matrix[m, m] Psi_star; // original was: L_Psi);
    matrix[n, m] Pi_t;
    matrix[m, p] L_Yt;
    matrix[p, m] L_Y_A[Ng];
    matrix[n, q] L_Xt;
    matrix[n, m] cov_eta_xi;
    matrix[q, m] cov_x_eta;
    matrix[n, p] cov_y_xi;
    matrix[q, p] cov_y_x;
    matrix[n, q] cov_x_xi;
    matrix[m, m] cov_eta;
    
    matrix[p + q, p + q] top_left;
      
    matrix[m + n, p + q] corner;
    
    matrix[m + n, m + n] bottom_right;
    
    matrix[p + q, p + q] precision;
    matrix[m + n, m + n] L;
    matrix[m + n, p + q] beta;
    vector[m + n] lvmean;
    vector[p + q] ovmean[Ng];

    int obsidx[p + q];
    int r1 = 1;
    int r2 = 1;
    int grpidx = 1;

    for (g in 1:Ng) {
      ovmean[g] = to_vector(Nu[g]);

      if (q > 0 && n > 0) {
	ovmean[g, (p + 1):(p + q)] += to_vector(Lambda_x[g] * Alpha[g, (m + 1):(m + n), 1]);
      }

      if (p > 0) {
	L_Y_A[g] = mdivide_right(L_Y[g], I - Bet[g]);
	if (m > 0) {
	  ovmean[g, 1:p] += to_vector(L_Y_A[g] * Alpha[g, 1:m, 1]);
	}
	if (n > 0) {
	  ovmean[g, 1:p] += to_vector(L_Y_A[g] * Gam[g] * Alpha[g, (m + 1):(m + n), 1]);
	}
      }
    }
    
    for (mm in 1:Np) {
      grpidx = grpnum[mm];

      A = mdivide_left_tri_low(I - Bet[grpidx], I); // = (I - B)^{-1}
      total_eta_eta = A - I;
      indirect_eta_eta = total_eta_eta - Bet[grpidx];
      total_eta_y = L_Y[grpidx] * A;
      indirect_eta_y = total_eta_y - L_Y[grpidx];

      if (n > 0) {
	total_xi_eta = A * Gam[grpidx];
	indirect_xi_eta = total_xi_eta - Gam[grpidx];
	total_xi_y = total_eta_y * Gam[grpidx];
      }

      Psi_star = multiply_lower_tri_self_transpose(A * Psi[grpidx]); // original was: L_Psi);
      Pi_t = transpose(total_xi_eta);
      L_Yt = transpose(L_Y[grpidx]);
      L_Xt = transpose(L_X[grpidx]);

      if (n > 0) {
        cov_eta_xi = PHI[grpidx] * Pi_t;
        cov_x_eta = L_X[grpidx] * cov_eta_xi;
        cov_y_xi = cov_eta_xi * L_Yt;
        cov_y_x = L_X[grpidx] * cov_y_xi;
        cov_x_xi = PHI[grpidx] * L_Xt;
        cov_eta = quad_form_sym(PHI[grpidx], Pi_t) + Psi_star;

        top_left = append_row(
          append_col(quad_form_sym(cov_eta, L_Yt) + Theta[grpidx], 
                     transpose(cov_y_x)), 
          append_col(cov_y_x, quad_form_sym(PHI[grpidx], L_Xt) + Theta_x[grpidx]) );
      
        corner = transpose(append_col(
          append_row(cov_y_xi, cov_eta * L_Yt),
          append_row(cov_x_xi, transpose(cov_x_eta))));
    
        bottom_right = append_row(
          append_col(cov_eta, transpose(cov_eta_xi)), append_col(cov_eta_xi, PHI[grpidx]) );
      } else {
        cov_eta = Psi_star;
        top_left = quad_form_sym(cov_eta, L_Yt) + Theta[grpidx];
      
        corner = cov_eta * L_Yt;
    
        bottom_right = cov_eta;
      }

      // FIXME?? what if obsidx also extends to x variables?
      obsidx = Obsvar[mm, ];
      precision[1:Nobs[mm], 1:Nobs[mm]] = inverse_spd(top_left[obsidx[1:Nobs[mm]], obsidx[1:Nobs[mm]]]);
      L = cholesky_decompose(bottom_right - quad_form(precision[1:Nobs[mm], 1:Nobs[mm]], transpose(corner[, obsidx[1:Nobs[mm]]])));
      beta[, 1:Nobs[mm]] = corner[, obsidx[1:Nobs[mm]]] * precision[1:Nobs[mm], 1:Nobs[mm]];

      r1 = startrow[mm];
      r2 = endrow[mm];

      for (idx in r1:r2){
	lvmean = Alpha[grpidx, , 1] + beta[, 1:Nobs[mm]] * (YX[idx, 1:Nobs[mm]] - ovmean[grpidx, obsidx[1:Nobs[mm]]]);
	eta[idx,] = transpose(multi_normal_cholesky_rng(lvmean, L));
      }
    }
  }
} // end a with a completely blank line (not even whitespace)
