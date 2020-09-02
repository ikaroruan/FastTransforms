void X(destroy_tdc_eigen)(X(tdc_eigen) * F) {
    if (F->n < TDC_EIGEN_BLOCKSIZE) {
        free(F->V);
        free(F->lambda);
    }
    else {
        X(destroy_symmetric_dpr1_eigen)(F->F0);
        X(destroy_tdc_eigen)(F->F1);
        X(destroy_tdc_eigen)(F->F2);
        free(F->z);
    }
    free(F);
}

void X(destroy_tdc_eigen_FMM)(X(tdc_eigen_FMM) * F) {
    if (F->n < TDC_EIGEN_BLOCKSIZE) {
        free(F->V);
        free(F->lambda);
    }
    else {
        X(destroy_symmetric_dpr1_eigen_FMM)(F->F0);
        X(destroy_tdc_eigen_FMM)(F->F1);
        X(destroy_tdc_eigen_FMM)(F->F2);
        free(F->z);
    }
    free(F);
}

size_t X(summary_size_tdc_eigen)(X(tdc_eigen) * F) {
    size_t S = 0;
    if (F->n < TDC_EIGEN_BLOCKSIZE)
        S += sizeof(FLT)*F->n*(F->n+1);
    else {
        S += X(summary_size_symmetric_dpr1_eigen)(F->F0);
        S += X(summary_size_tdc_eigen)(F->F1);
        S += X(summary_size_tdc_eigen)(F->F2);
    }
    return S;
}

size_t X(summary_size_tdc_eigen_FMM)(X(tdc_eigen_FMM) * F) {
    size_t S = 0;
    if (F->n < TDC_EIGEN_BLOCKSIZE)
        S += sizeof(FLT)*F->n*(F->n+1);
    else {
        S += X(summary_size_symmetric_dpr1_eigen_FMM)(F->F0);
        S += X(summary_size_tdc_eigen_FMM)(F->F1);
        S += X(summary_size_tdc_eigen_FMM)(F->F2);
    }
    return S;
}



void Y(symmetric_tridiagonal_printmat)(char * MAT, char * FMT, X(symmetric_tridiagonal) * A) {
    int n = A->n;
    FLT * a = A->a;
    FLT * b = A->b;

    printf("%s = \n", MAT);
    if (n == 1) {
        if (a[0] < 0) {printf("[");}
        else {printf("[ ");}
        printf(FMT, a[0]);
    }
    else if (n == 2) {
        if (a[0] < 0) {printf("[");}
        else {printf("[ ");}
        printf(FMT, a[0]);
        if (b[0] < 0) {printf("  ");}
        else {printf("   ");}
        printf(FMT, b[0]);
        printf("\n");
        if (b[0] < 0) {printf(" ");}
        else {printf("  ");}
        printf(FMT, b[0]);
        if (a[1] < 0) {printf("  ");}
        else {printf("   ");}
        printf(FMT, a[1]);
    }
    else if (n > 2) {
        // First row
        if (a[0] < 0) {printf("[");}
        else {printf("[ ");}
        printf(FMT, a[0]);
        if (b[0] < 0) {printf("  ");}
        else {printf("   ");}
        printf(FMT, b[0]);
        for (int j = 2; j < n; j++) {
            printf("   ");
            printf(FMT, 0.0);
        }

        // Second row
        printf("\n");
        if (b[0] < 0) {printf(" ");}
        else {printf("  ");}
        printf(FMT, b[0]);
        if (a[1] < 0) {printf("  ");}
        else {printf("   ");}
        printf(FMT, a[1]);
        if (b[1] < 0) {printf("  ");}
        else {printf("   ");}
        printf(FMT, b[1]);
        for (int j = 3; j < n; j++) {
            printf("   ");
            printf(FMT, 0.0);
        }

        // Interior rows
        for (int i = 2; i < n-1; i++) {
            printf("\n");
            printf("  ");
            printf(FMT, 0.0);
            for (int j = 1; j < i-1; j++) {
                printf("   ");
                printf(FMT, 0.0);
            }
            if (b[i-1] < 0) {printf("  ");}
            else {printf("   ");}
            printf(FMT, b[i-1]);
            if (a[i] < 0) {printf("  ");}
            else {printf("   ");}
            printf(FMT, a[i]);
            if (b[i] < 0) {printf("  ");}
            else {printf("   ");}
            printf(FMT, b[i]);
            for (int j = i+2; j < n; j++) {
                printf("   ");
                printf(FMT, 0.0);
            }
        }

        // Last row
        printf("\n");
        printf("  ");
        printf(FMT, 0.0);
        for (int j = 1; j < n-2; j++) {
            printf("   ");
            printf(FMT, 0.0);
        }
        if (b[n-2] < 0) {printf("  ");}
        else {printf("   ");}
        printf(FMT, b[n-2]);
        if (a[n-1] < 0) {printf("  ");}
        else {printf("   ");}
        printf(FMT, a[n-1]);
    }
    printf("]\n");
}

#define A(i,j) A[(i)+n*(j)]
#define B(i,j) B[(i)+n*(j)]

void X(printmat)(char * MAT, char * FMT, FLT * A, int n, int m) {
    printf("%s = \n", MAT);
    if (n > 0 && m > 0) {
        if (signbit(A(0,0))) {printf("[");}
        else {printf("[ ");}
        printf(FMT, A(0,0));
        for (int j = 1; j < m; j++) {
            if (signbit(A(0,j))) {printf("  ");}
            else {printf("   ");}
            printf(FMT, A(0,j));
        }
        for (int i = 1; i < n-1; i++) {
            printf("\n");
            if (signbit(A(i,0))) {printf(" ");}
            else {printf("  ");}
            printf(FMT, A(i,0));
            for (int j = 1; j < m; j++) {
                if (signbit(A(i,j))) {printf("  ");}
                else {printf("   ");}
                printf(FMT, A(i,j));
            }
        }
        if (n > 1) {
            printf("\n");
            if (signbit(A(n-1,0))) {printf(" ");}
            else {printf("  ");}
            printf(FMT, A(n-1,0));
            for (int j = 1; j < m; j++) {
                if (signbit(A(n-1,j))) {printf("  ");}
                else {printf("   ");}
                printf(FMT, A(n-1,j));
            }
        }
        printf("]\n");
    }
}

X(tdc_eigen) * X(tdc_eig)(X(symmetric_tridiagonal) * A) {
	int n = A->n;
	X(tdc_eigen) * F = malloc(sizeof(X(tdc_eigen)));
	
	if(n < TDC_EIGEN_BLOCKSIZE){
		FLT * V = calloc(n*n, sizeof(FLT));
		for(int i = 0; i < n; i++)
			V[i + i*n] = 1;
		FLT * lambda = calloc(n, sizeof(FLT));

		X(symmetric_tridiagonal_eig)(A, V, lambda);
		F->V = V;
		F->lambda = lambda;
		F->n = n;
	}
	else{
		int s = n >> 1;
		int sign = -1;
		FLT rho = -sign*Y(fabs)(A->b[s-1]);

		// Perturbed blocks.
		X(symmetric_tridiagonal) * A1t = malloc(sizeof(X(symmetric_tridiagonal)));
		FLT * A1ta = malloc(s*sizeof(FLT));
		FLT * A1tb = malloc((s-1)*sizeof(FLT));
		for(int i = 0; i < s-1; i++){
			A1ta[i] = A->a[i];
			A1tb[i] = A->b[i];
		}
		A1ta[s-1] = A->a[s-1] + sign*Y(fabs)(rho);
		A1t->a = A1ta;
		A1t->b = A1tb;
		A1t->n = s;

		X(symmetric_tridiagonal) * A2t = malloc(sizeof(X(symmetric_tridiagonal)));
		FLT * A2ta = malloc((n-s)*sizeof(FLT));
		FLT * A2tb = malloc((n-s-1)*sizeof(FLT));
		for(int i = 0; i < n-s-1; i++){
			A2ta[i+1] = A->a[s+i+1];
			A2tb[i] = A->b[s+i];
		}
		A2ta[0] = A->a[s] + sign*Y(fabs)(rho);
		A2t->a = A2ta;
		A2t->b = A2tb;
		A2t->n = n-s;

		F->F1 = X(tdc_eig)(A1t);
		F->F2 = X(tdc_eig)(A2t);

		FLT * y = calloc(n, sizeof(FLT));
		y[s-1] = sign; y[s] = 1;
		FLT * z = calloc(n, sizeof(FLT));
		X(tdmv)('T', 1, F->F1, y, 0, z);
		X(tdmv)('T', 1, F->F2, y+s, 0, z+s);

		FLT * d = malloc(n*sizeof(FLT));
		for(int i = 0; i < s; i++){
			if(F->F1->n < TDC_EIGEN_BLOCKSIZE)
				d[i] = F->F1->lambda[i];
			else
				d[i] = F->F1->F0->lambda[i];
		}
		for(int i = 0; i < n-s; i++){
			if(F->F2->n < TDC_EIGEN_BLOCKSIZE)
				d[i+s] = F->F2->lambda[i];
			else
				d[i+s] = F->F2->F0->lambda[i];
		}

		X(symmetric_dpr1) * Ar = malloc(sizeof(X(symmetric_dpr1)));
		Ar->d = d;
		Ar->z = z;
		Ar->rho = rho;
		Ar->n = n;

		F->F0 = X(symmetric_dpr1_eig)(Ar);
		F->n = n;
		F->z = y;

		X(destroy_symmetric_tridiagonal)(A1t);
		X(destroy_symmetric_tridiagonal)(A2t);
		X(destroy_symmetric_dpr1)(Ar);
	}

	return F;
}

X(tdc_eigen_FMM) * X(tdc_eig_FMM)(X(symmetric_tridiagonal) * A) {
	int n = A->n;
	X(tdc_eigen_FMM) * F = malloc(sizeof(X(tdc_eigen_FMM)));
	
	if(n < TDC_EIGEN_BLOCKSIZE){
		FLT * V = calloc(n*n, sizeof(FLT));
		for(int i = 0; i < n; i++)
			V[i + i*n] = 1;
		FLT * lambda = calloc(n, sizeof(FLT));

		X(symmetric_tridiagonal_eig)(A, V, lambda);
		F->V = V;
		F->lambda = lambda;
		F->n = n;
	}
	else{
		int s = n >> 1;
		int sign = -1;
		FLT rho = -sign*Y(fabs)(A->b[s-1]);

		// Perturbed blocks.
		X(symmetric_tridiagonal) * A1t = malloc(sizeof(X(symmetric_tridiagonal)));
		FLT * A1ta = malloc(s*sizeof(FLT));
		FLT * A1tb = malloc((s-1)*sizeof(FLT));
		for(int i = 0; i < s-1; i++){
			A1ta[i] = A->a[i];
			A1tb[i] = A->b[i];
		}
		A1ta[s-1] = A->a[s-1] + sign*Y(fabs)(rho);
		A1t->a = A1ta;
		A1t->b = A1tb;
		A1t->n = s;

		X(symmetric_tridiagonal) * A2t = malloc(sizeof(X(symmetric_tridiagonal)));
		FLT * A2ta = malloc((n-s)*sizeof(FLT));
		FLT * A2tb = malloc((n-s-1)*sizeof(FLT));
		for(int i = 0; i < n-s-1; i++){
			A2ta[i+1] = A->a[s+i+1];
			A2tb[i] = A->b[s+i];
		}
		A2ta[0] = A->a[s] + sign*Y(fabs)(rho);
		A2t->a = A2ta;
		A2t->b = A2tb;
		A2t->n = n-s;

		F->F1 = X(tdc_eig_FMM)(A1t);
		F->F2 = X(tdc_eig_FMM)(A2t);

		FLT * y = calloc(n, sizeof(FLT));
		y[s-1] = sign; y[s] = 1;
		FLT * z = calloc(n, sizeof(FLT));
		X(tfmv)('T', 1, F->F1, y, 0, z);
		X(tfmv)('T', 1, F->F2, y+s, 0, z+s);

		FLT * d = malloc(n*sizeof(FLT));
		for(int i = 0; i < s; i++){
			if(F->F1->n < TDC_EIGEN_BLOCKSIZE)
				d[i] = F->F1->lambda[i];
			else
				d[i] = F->F1->F0->lambda[i];
		}
		for(int i = 0; i < n-s; i++){
			if(F->F2->n < TDC_EIGEN_BLOCKSIZE)
				d[i+s] = F->F2->lambda[i];
			else
				d[i+s] = F->F2->F0->lambda[i];
		}

		X(symmetric_dpr1) * Ar = malloc(sizeof(X(symmetric_dpr1)));
		Ar->d = d;
		Ar->z = z;
		Ar->rho = rho;
		Ar->n = n;

		F->F0 = X(symmetric_dpr1_eig_FMM)(Ar);
		F->n = n;
		F->z = y;

		X(destroy_symmetric_tridiagonal)(A1t);
		X(destroy_symmetric_tridiagonal)(A2t);
		X(destroy_symmetric_dpr1)(Ar);
	}

	return F;
}

X(tdc_eigen) * X(sdtdc_eig)(X(symmetric_tridiagonal) * T, X(symmetric_tridiagonal) * S) {
    int n = T->n;
    X(tdc_eigen) * F = malloc(sizeof(X(tdc_eigen)));
    if (n < TDC_EIGEN_BLOCKSIZE) {
        FLT * V = calloc(n*n, sizeof(FLT));
        for (int i = 0; i < n; i++)
            V[i+i*n] = 1;
        FLT * lambda = calloc(n, sizeof(FLT));
        X(symmetric_definite_tridiagonal_eig)(T, S, V, lambda);
        F->V = V;
        F->lambda = lambda;
        F->n = n;
    }
    else {
        int s = n>>1, sign = -1;

        FLT * y = calloc(n, sizeof(FLT));
        y[s-1] = sign;
        y[s] = 1;

        FLT rho = -sign*Y(fabs)(T->b[s-1]);
        FLT sigma = -sign*Y(fabs)(S->b[s-1]);

        X(symmetric_tridiagonal) * T1 = malloc(sizeof(X(symmetric_tridiagonal)));
        FLT * at1 = malloc(s*sizeof(FLT));
        FLT * bt1 = malloc((s-1)*sizeof(FLT));
        for (int i = 0; i < s-1; i++) {
            at1[i] = T->a[i];
            bt1[i] = T->b[i];
        }
        at1[s-1] = T->a[s-1] + sign*Y(fabs)(rho);
        T1->a = at1;
        T1->b = bt1;
        T1->n = s;

        X(symmetric_tridiagonal) * S1 = malloc(sizeof(X(symmetric_tridiagonal)));
        FLT * as1 = malloc(s*sizeof(FLT));
        FLT * bs1 = malloc((s-1)*sizeof(FLT));
        for (int i = 0; i < s-1; i++) {
            as1[i] = S->a[i];
            bs1[i] = S->b[i];
        }
        as1[s-1] = S->a[s-1] + sign*Y(fabs)(sigma);
        S1->a = as1;
        S1->b = bs1;
        S1->n = s;

        X(symmetric_tridiagonal) * T2 = malloc(sizeof(X(symmetric_tridiagonal)));
        FLT * at2 = malloc((n-s)*sizeof(FLT));
        FLT * bt2 = malloc((n-s-1)*sizeof(FLT));
        for (int i = s; i < n-1; i++) {
            at2[i-s+1] = T->a[i+1];
            bt2[i-s] = T->b[i];
        }
        at2[0] = T->a[s] + sign*Y(fabs)(rho);
        T2->a = at2;
        T2->b = bt2;
        T2->n = n-s;

        X(symmetric_tridiagonal) * S2 = malloc(sizeof(X(symmetric_tridiagonal)));
        FLT * as2 = malloc((n-s)*sizeof(FLT));
        FLT * bs2 = malloc((n-s-1)*sizeof(FLT));
        for (int i = s; i < n-1; i++) {
            as2[i-s+1] = S->a[i+1];
            bs2[i-s] = S->b[i];
        }
        as2[0] = S->a[s] + sign*Y(fabs)(sigma);
        S2->a = as2;
        S2->b = bs2;
        S2->n = n-s;

        F->F1 = X(sdtdc_eig)(T1, S1);
        F->F2 = X(sdtdc_eig)(T2, S2);

        FLT * z = calloc(n, sizeof(FLT));
        X(tdmv)('T', 1, F->F1, y, 0, z);
        X(tdmv)('T', 1, F->F2, y+s, 0, z+s);

        FLT * d = malloc(n*sizeof(FLT));
        for (int i = 0; i < s; i++) {
            if (F->F1->n < TDC_EIGEN_BLOCKSIZE)
                d[i] = F->F1->lambda[i];
            else
                d[i] = F->F1->F0->lambda[i];
        }
        for (int i = 0; i < n-s; i++) {
            if (F->F2->n < TDC_EIGEN_BLOCKSIZE)
                d[i+s] = F->F2->lambda[i];
            else
                d[i+s] = F->F2->F0->lambda[i];
        }

        X(symmetric_dpr1) * A = malloc(sizeof(X(symmetric_dpr1)));
        A->d = d;
        A->z = z;
        A->rho = rho;
        A->n = n;
        FLT * zb = malloc(n*sizeof(FLT));
        for (int i = 0; i < n; i++)
            zb[i] = z[i];
        X(symmetric_idpr1) * B = malloc(sizeof(X(symmetric_idpr1)));
        B->z = zb;
        B->sigma = sigma;
        B->n = n;

        F->F0 = X(symmetric_definite_dpr1_eig)(A, B);
        F->n = n;

        X(destroy_symmetric_tridiagonal)(T1);
        X(destroy_symmetric_tridiagonal)(T2);
        X(destroy_symmetric_tridiagonal)(S1);
        X(destroy_symmetric_tridiagonal)(S2);
        X(destroy_symmetric_dpr1)(A);
        X(destroy_symmetric_idpr1)(B);
        F->z = y;
    }
    return F;
}

X(tdc_eigen_FMM) * X(sdtdc_eig_FMM)(X(symmetric_tridiagonal) * T, X(symmetric_tridiagonal) * S) {
    int n = T->n;
    X(tdc_eigen_FMM) * F = malloc(sizeof(X(tdc_eigen_FMM)));
    if (n < TDC_EIGEN_BLOCKSIZE) {
        FLT * V = calloc(n*n, sizeof(FLT));
        for (int i = 0; i < n; i++)
            V[i+i*n] = 1;
        FLT * lambda = calloc(n, sizeof(FLT));
        X(symmetric_definite_tridiagonal_eig)(T, S, V, lambda);
        F->V = V;
        F->lambda = lambda;
        F->n = n;
    }
    else {
        int s = n>>1, sign = -1;

        FLT * y = calloc(n, sizeof(FLT));
        y[s-1] = sign;
        y[s] = 1;

        FLT rho = -sign*Y(fabs)(T->b[s-1]);
        FLT sigma = -sign*Y(fabs)(S->b[s-1]);

        X(symmetric_tridiagonal) * T1 = malloc(sizeof(X(symmetric_tridiagonal)));
        FLT * at1 = malloc(s*sizeof(FLT));
        FLT * bt1 = malloc((s-1)*sizeof(FLT));
        for (int i = 0; i < s-1; i++) {
            at1[i] = T->a[i];
            bt1[i] = T->b[i];
        }
        at1[s-1] = T->a[s-1] + sign*Y(fabs)(rho);
        T1->a = at1;
        T1->b = bt1;
        T1->n = s;

        X(symmetric_tridiagonal) * S1 = malloc(sizeof(X(symmetric_tridiagonal)));
        FLT * as1 = malloc(s*sizeof(FLT));
        FLT * bs1 = malloc((s-1)*sizeof(FLT));
        for (int i = 0; i < s-1; i++) {
            as1[i] = S->a[i];
            bs1[i] = S->b[i];
        }
        as1[s-1] = S->a[s-1] + sign*Y(fabs)(sigma);
        S1->a = as1;
        S1->b = bs1;
        S1->n = s;

        X(symmetric_tridiagonal) * T2 = malloc(sizeof(X(symmetric_tridiagonal)));
        FLT * at2 = malloc((n-s)*sizeof(FLT));
        FLT * bt2 = malloc((n-s-1)*sizeof(FLT));
        for (int i = s; i < n-1; i++) {
            at2[i-s+1] = T->a[i+1];
            bt2[i-s] = T->b[i];
        }
        at2[0] = T->a[s] + sign*Y(fabs)(rho);
        T2->a = at2;
        T2->b = bt2;
        T2->n = n-s;

        X(symmetric_tridiagonal) * S2 = malloc(sizeof(X(symmetric_tridiagonal)));
        FLT * as2 = malloc((n-s)*sizeof(FLT));
        FLT * bs2 = malloc((n-s-1)*sizeof(FLT));
        for (int i = s; i < n-1; i++) {
            as2[i-s+1] = S->a[i+1];
            bs2[i-s] = S->b[i];
        }
        as2[0] = S->a[s] + sign*Y(fabs)(sigma);
        S2->a = as2;
        S2->b = bs2;
        S2->n = n-s;

        F->F1 = X(sdtdc_eig_FMM)(T1, S1);
        F->F2 = X(sdtdc_eig_FMM)(T2, S2);

        FLT * z = calloc(n, sizeof(FLT));
        X(tfmv)('T', 1, F->F1, y, 0, z);
        X(tfmv)('T', 1, F->F2, y+s, 0, z+s);

        FLT * d = malloc(n*sizeof(FLT));
        for (int i = 0; i < s; i++) {
            if (F->F1->n < TDC_EIGEN_BLOCKSIZE)
                d[i] = F->F1->lambda[i];
            else
                d[i] = F->F1->F0->lambda[i];
        }
        for (int i = 0; i < n-s; i++) {
            if (F->F2->n < TDC_EIGEN_BLOCKSIZE)
                d[i+s] = F->F2->lambda[i];
            else
                d[i+s] = F->F2->F0->lambda[i];
        }

        X(symmetric_dpr1) * A = malloc(sizeof(X(symmetric_dpr1)));
        A->d = d;
        A->z = z;
        A->rho = rho;
        A->n = n;
        FLT * zb = malloc(n*sizeof(FLT));
        for (int i = 0; i < n; i++)
            zb[i] = z[i];
        X(symmetric_idpr1) * B = malloc(sizeof(X(symmetric_idpr1)));
        B->z = zb;
        B->sigma = sigma;
        B->n = n;

        F->F0 = X(symmetric_definite_dpr1_eig_FMM)(A, B);
        F->n = n;

        X(destroy_symmetric_tridiagonal)(T1);
        X(destroy_symmetric_tridiagonal)(T2);
        X(destroy_symmetric_tridiagonal)(S1);
        X(destroy_symmetric_tridiagonal)(S2);
        X(destroy_symmetric_dpr1)(A);
        X(destroy_symmetric_idpr1)(B);
        F->z = y;
    }
    return F;
}

void X(tdmv)(char TRANS, FLT alpha, X(tdc_eigen) * F, FLT * x, FLT beta, FLT * y) {
    int n = F->n;
    if (n < TDC_EIGEN_BLOCKSIZE)
        X(gemv)(TRANS, n, n, alpha, F->V, n, x, beta, y);
    else {
        int s = n>>1;
        FLT * z = F->z;
        if (TRANS == 'N') {
            X(dvmv)(TRANS, 1, F->F0, x, 0, z);
            X(tdmv)(TRANS, alpha, F->F1, z, beta, y);
            X(tdmv)(TRANS, alpha, F->F2, z+s, beta, y+s);
        }
        else if (TRANS == 'T') {
            X(tdmv)(TRANS, 1, F->F1, x, 0, z);
            X(tdmv)(TRANS, 1, F->F2, x+s, 0, z+s);
            X(dvmv)(TRANS, alpha, F->F0, z, beta, y);
        }
    }
}

void X(tfmv)(char TRANS, FLT alpha, X(tdc_eigen_FMM) * F, FLT * x, FLT beta, FLT * y) {
    int n = F->n;
    if (n < TDC_EIGEN_BLOCKSIZE)
        X(gemv)(TRANS, n, n, alpha, F->V, n, x, beta, y);
    else {
        int s = n>>1;
        FLT * z = F->z;
        if (TRANS == 'N') {
            X(dfmv)(TRANS, 1, F->F0, x, 0, z);
            X(tfmv)(TRANS, alpha, F->F1, z, beta, y);
            X(tfmv)(TRANS, alpha, F->F2, z+s, beta, y+s);
        }
        else if (TRANS == 'T') {
            X(tfmv)(TRANS, 1, F->F1, x, 0, z);
            X(tfmv)(TRANS, 1, F->F2, x+s, 0, z+s);
            X(dfmv)(TRANS, alpha, F->F0, z, beta, y);
        }
    }
}
