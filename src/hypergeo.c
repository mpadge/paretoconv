#include <math.h>
#include <float.h>
#include <R.h>
#include <Rinternals.h>

#define ARRAY_LENGTH 777
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/* ===== Exponential Form Operations ===== */

/* EMULT: Multiply two numbers in exponential form (base, exponent) */
void emult(double n1, double e1, double n2, double e2,
           double *nf, double *ef) {
    *nf = n1 * n2;
    *ef = e1 + e2;
    if (fabs(*nf) >= 10.0) {
        *nf = *nf / 10.0;
        *ef = *ef + 1.0;
    }
}

/* EDIV: Divide two numbers in exponential form */
void ediv(double n1, double e1, double n2, double e2,
          double *nf, double *ef) {
    *nf = n1 / n2;
    *ef = e1 - e2;
    if ((fabs(*nf) < 1.0) && (*nf != 0.0)) {
        *nf = *nf * 10.0;
        *ef = *ef - 1.0;
    }
}

/* EADD: Add two numbers in exponential form */
void eadd(double n1, double e1, double n2, double e2,
          double *nf, double *ef) {
    double ediff = e1 - e2;

    if (ediff > 36.0) {
        *nf = n1;
        *ef = e1;
    } else if (ediff < -36.0) {
        *nf = n2;
        *ef = e2;
    } else {
        *nf = n1 * pow(10.0, ediff) + n2;
        *ef = e2;
        /* Normalize */
        while (fabs(*nf) >= 10.0) {
            *nf = *nf / 10.0;
            *ef = *ef + 1.0;
        }
        while ((fabs(*nf) < 1.0) && (*nf != 0.0)) {
            *nf = *nf * 10.0;
            *ef = *ef - 1.0;
        }
    }
}

/* ESUB: Subtract two numbers in exponential form */
void esub(double n1, double e1, double n2, double e2,
          double *nf, double *ef) {
    eadd(n1, e1, -n2, e2, nf, ef);
}

/* ===== Array Operations ===== */

/* ARADD: Add two arbitrary-precision arrays
   Arrays use indexing: index 0 for sign, 1..L for mantissa, L+1 for exponent */
void aradd(double a[], double b[], double c[], int l, double rmax) {
    double z[ARRAY_LENGTH + 3];
    int i, ediff;

    for (i = -1; i <= l + 1; i++) {
        z[i + 1] = 0.0;
    }

    ediff = (int)floor(a[l + 2] - b[l + 2] + 0.5);

    if (fabs(a[0]) < 0.5 || ediff <= -l) {
        for (i = -1; i <= l + 1; i++) {
            c[i + 1] = b[i + 1];
        }
        return;
    }
    if (fabs(b[0]) < 0.5 || ediff >= l) {
        for (i = -1; i <= l + 1; i++) {
            c[i + 1] = a[i + 1];
        }
        return;
    }

    z[0] = a[0];

    if (fabs(a[0] - b[0]) < 0.5) {
        z[l + 2] = a[l + 2];
        for (i = l; i >= 1; i--) {
            z[i + 1] = a[i + 1] + b[i + 1] + z[i + 1];
            if (z[i + 1] >= rmax) {
                z[i + 1] = z[i + 1] - rmax;
                z[i] = 1.0;
            }
        }
    } else {
        if (ediff > 0) {
            z[l + 2] = a[l + 2];
            for (i = l; i >= 1; i--) {
                z[i + 1] = a[i + 1] - b[i + 1] + z[i + 1];
                if (z[i + 1] < 0.0) {
                    z[i + 1] = z[i + 1] + rmax;
                    z[i] = -1.0;
                }
            }
        } else {
            z[l + 2] = b[l + 2];
            z[0] = b[0];
            for (i = l; i >= 1; i--) {
                z[i + 1] = b[i + 1] - a[i + 1] + z[i + 1];
                if (z[i + 1] < 0.0) {
                    z[i + 1] = z[i + 1] + rmax;
                    z[i] = -1.0;
                }
            }
        }
    }

    for (i = -1; i <= l + 1; i++) {
        c[i + 1] = z[i + 1];
    }
}

/* ARSUB: Subtract two arrays */
void arsub(double a[], double b[], double c[], int l, double rmax) {
    double b2[ARRAY_LENGTH + 3];
    int i;

    for (i = -1; i <= l + 1; i++) {
        b2[i + 1] = b[i + 1];
    }
    b2[0] = -b2[0];
    aradd(a, b2, c, l, rmax);
}

/* ARMULT: Multiply array by scalar */
void armult(double a[], double b, double c[], int l, double rmax) {
    double z[ARRAY_LENGTH + 3];
    double carry;
    int i;

    z[0] = copysign(1.0, b) * a[0];
    double b2 = fabs(b);
    z[l + 2] = a[l + 2];

    for (i = 0; i <= l; i++) {
        z[i + 1] = 0.0;
    }

    if (b2 <= 1.0e-10 || a[1] <= 1.0e-10) {
        z[0] = 1.0;
        z[l + 2] = 0.0;
    } else {
        for (i = l; i >= 1; i--) {
            z[i + 1] = a[i + 1] * b2 + z[i + 1];
            if (z[i + 1] >= rmax) {
                carry = floor(z[i + 1] / rmax);
                z[i + 1] = z[i + 1] - carry * rmax;
                z[i] = carry;
            }
        }
        if (z[1] >= 0.5) {
            for (i = l; i >= 1; i--) {
                z[i + 1] = z[i];
            }
            z[l + 2] = z[l + 2] + 1.0;
            z[1] = 0.0;
        }
    }

    for (i = -1; i <= l + 1; i++) {
        c[i + 1] = z[i + 1];
    }
}

/* ===== Main Hypergeometric Function ===== */

/* CHGF: Confluent hypergeometric function 1F1(a,b,z) for real arguments */
double chgf(double a, double b, double z, int l) {
    int bit = DBL_MANT_DIG;
    double rmax = pow(2.0, (double)bit / 2.0);
    double sigfig = pow(2.0, (double)bit / 4.0);

    double sum[ARRAY_LENGTH + 3], num[ARRAY_LENGTH + 3], denom[ARRAY_LENGTH + 3];
    double q1[ARRAY_LENGTH + 3], q2[ARRAY_LENGTH + 3];

    /* Extract components with high precision */
    double ar, ar2;
    ar2 = a * sigfig;
    ar = floor(ar2);
    ar2 = nearbyint((ar2 - ar) * rmax);

    double cr, cr2;
    cr2 = b * sigfig;
    cr = floor(cr2);
    cr2 = nearbyint((cr2 - cr) * rmax);

    double xr, xr2;
    xr2 = z * sigfig;
    xr = floor(xr2);
    xr2 = nearbyint((xr2 - xr) * rmax);

    /* Initialize arrays */
    int i;
    for (i = -1; i <= l + 1; i++) {
        sum[i + 1] = 0.0;
        num[i + 1] = 0.0;
        denom[i + 1] = 0.0;
    }

    sum[0] = 1.0;
    num[0] = 1.0;
    denom[0] = 1.0;

    sum[2] = 1.0;
    num[2] = 1.0;
    denom[2] = 1.0;

    double cnt = sigfig;

    /* Kummer series summation */
    while (1) {
        if (sum[1] < 0.5) {
            if (fabs(sum[l + 2]) > 2.0) {
                if (cr > 0.0) {
                    double term_abs = fabs(ar + cr) * fabs(xr) / (fabs(cr + cr2) * cnt);
                    if (term_abs <= 1.0) {
                        break;
                    }
                }
            }
        }

        /* sum *= (cr + cr2) * xr / (cr * cnt) */
        armult(sum, xr, q1, l, rmax);
        armult(q1, cr, q2, l, rmax);
        armult(q2, cr2, sum, l, rmax);
        armult(sum, 1.0 / cnt, sum, l, rmax);

        /* denom *= (cr + cr2) * cnt */
        armult(denom, cr, q1, l, rmax);
        armult(q1, cr2, denom, l, rmax);
        armult(denom, cnt, denom, l, rmax);

        /* num *= (ar + cnt) * xr */
        armult(num, ar + cnt, q1, l, rmax);
        armult(q1, xr, num, l, rmax);

        /* sum += num / denom (approximate) */
        if (num[1] > 0.5 && denom[1] > 0.5) {
            aradd(sum, num, sum, l, rmax);
        }

        cnt += sigfig;
        ar += sigfig;
        cr += sigfig;
    }

    /* Final division: extract result */
    if (denom[1] < 0.5) {
        return copysign(pow(10.0, sum[l + 2]), sum[0]);
    }

    double n1, e1, n2, e2, n3, e3;
    emult(sum[1], sum[l + 2], sum[1], sum[l + 2], &n1, &e1);
    emult(denom[1], denom[l + 2], denom[1], denom[l + 2], &n2, &e2);
    eadd(n1, e1, n2, e2, &n3, &e3);

    ediv(sum[1], sum[l + 2], denom[1], denom[l + 2], &n1, &e1);
    return copysign(n1 * pow(10.0, e1), sum[0]);
}

/* CONHYP: Determine appropriate parameters and call CHGF */
double conhyp(double a, double b, double z) {
    int b_val = DBL_MANT_DIG;

    double nterm = 0;
    double fx = 0;
    double term1 = 0;

    while (1) {
        nterm += 1.0;
        double term_num = (a + nterm - 1.0) * z;
        double term_den = (b + nterm - 1.0) * nterm;
        double term2 = fabs(term_num / term_den);

        if (term2 == 0.0) {
            break;
        }

        if (term2 < 1.0) {
            if ((a + nterm - 1.0) > 1.0) {
                if ((b + nterm - 1.0) > 1.0) {
                    if ((term2 - term1) < 0.0) {
                        break;
                    }
                }
            }
        }

        fx += log(fmax(term2, 1.0e-10));
        term1 = term2;
    }

    double ang = (fabs(z) > 0) ? atan2(0.0, z) : 1.0;
    if (fabs(ang) < (M_PI * 0.5)) {
        ang = 1.0;
    } else {
        ang = sin(fabs(ang) - (M_PI * 0.5)) + 1.0;
    }

    double max = fx * 2.0 / (b_val * 6.93147181e-1);
    int i = (int)floor(max * ang) + 7;
    if (i < 5) i = 5;

    return chgf(a, b, z, i);
}

/* ===== R Interface ===== */

/* CHFM_CALL: R-callable wrapper for Kummer's confluent hypergeometric function
   Takes real arrays and returns real result vector */
SEXP chfm_call(SEXP z, SEXP a, SEXP b) {
    int n = LENGTH(z);

    double *z_vec = REAL(z);
    double a_val = asReal(a);
    double b_val = asReal(b);

    /* Create output vector */
    SEXP result = PROTECT(allocVector(REALSXP, n));
    double *result_vec = REAL(result);

    /* Compute hypergeometric function for each input */
    for (int i = 0; i < n; i++) {
        result_vec[i] = conhyp(a_val, b_val, z_vec[i]);
    }

    UNPROTECT(1);
    return result;
}
