/*
 * godunov.cc
 *
 * ������� ������ ������ ��� ������� ��������� ������� ��������.
 *
 * Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 3rd Edition. - Springer,
 * 2009. - P. 152 - 162.
 * 
 * (c) ����� �����, 2013
 *
 * ������: 26 ������� 2013 �.
 *
 */

#include "godunov.h"

#include "utils.h"

/* ������������ ��������� ������� �������� � �������� �� ���������� �������

   ���: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 3rd Edition. - Springer,
   2009. - P. 155. - Subroutine STARPU.

   params - ��������� � ����������� ��������������� ������������ (in)
   v_ncons_l[M] - ������ ����������� ���������� ����� �� ������� (in)
   v_ncons_r[M] - ������ ����������� ���������� ������ �� ������� (in)
   cl - �������� ����� ����� �� ������� (in)
   cr - �������� ����� ������ �� ������� (in)

   p_cont - �������� �� ���������� ������� (out)
   v_cont - �������� �� ���������� ������� (out) */
void calc_contact_pressure_velocity( struct Parameters *params, double v_ncons_l[M], double v_ncons_r[M],
                                     double cl, double cr, double *p_cont, double *v_cont ) {

    double vl, vr;      /* �������� ����� � ������ �� ������� */
    double p_old;       /* �������� �������� �� ���������� �������� */
    double fl, fr;      /* �������� ������� */
    double fld, frd;    /* �������� ����������� */
    int iter_num = 0;   /* ���������� ����������� �������� */
    double criteria;    /* ���������� ��� ����������� ���������� */
    double g;           /* ���������� �������� */

    /* �������� ����������� ��� �������� ������ ������ */
    vl = v_ncons_l[V];
    vr = v_ncons_r[V];
    g = params->g;

    if ( 2.0 * ( cl + cr ) / ( g - 1.0 ) <= vr - vl ) {
        /* ������ ������������� ������� */
        printf( "\ncalc_contact_pressure_velocity -> vacuum is generated\n" );
        exit( EXIT_FAILURE );
    }

    /* ������ ���������� ����������� ��� �������� */
    p_old = pressure_initial_guess( params, v_ncons_l, v_ncons_r, cl, cr );
    if ( p_old < 0.0 ) {
        printf( "\ncalc_contact_pressure_velocity -> initial pressure guess is negative " );
        exit( EXIT_FAILURE );
    }
    
    /* ������� ����������� ��������� ��� ���������� �������� �� ���������� ������� ������� �������-������� */
    do {
        calc_F_and_DF( params, p_old, v_ncons_l, cl, &fl, &fld );
        calc_F_and_DF( params, p_old, v_ncons_r, cr, &fr, &frd );
        *p_cont = p_old - ( fl + fr + vr - vl ) / ( fld + frd );
        criteria = 2.0 * fabs( ( *p_cont - p_old ) / ( *p_cont + p_old ) );
        iter_num++;
        if ( iter_num > MAX_ITER_NUM ) {
            printf( "\ncalc_contact_pressure_velocity -> number of iterations exceeds the maximum value.\n" );
            exit( EXIT_FAILURE );
        }
        if ( *p_cont < 0.0 ) {
            printf( "\ncalc_contact_pressure_velocity -> pressure is negative.\n" );
            exit( EXIT_FAILURE );            
        }
        p_old = *p_cont;
    } while ( criteria > EPS );

    /* �������� ����������� ������� */
    *v_cont = 0.5 * ( vl + vr + fr - fl );

}

/* ������ ������� F, ������������ �������� ���� �� ���������� �������, � �� ����������� �� �������� ����� DF

   ���: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 3rd Edition. - Springer,
   1999. - P. 158. - Subroutine PREFUN.
   + �������� �� ���������� ��������� ���������:

   params - ��������� � ����������� ��������������� ������������ (in)
   curr_press - �������� � ���������� �������� (in)
   v_ncons[M] - ������ ����������� ���������� (in)
   c - �������� ����� � ������������� ����� (in)

   F - �������� ������� (out)
   DF - �������� ����������� (out) */
void calc_F_and_DF( struct Parameters *params, double curr_press, double v_ncons[M], double c, double *F, double *DF ) {

    double r, v, p;                 /* ����������� ���������� */
    double g;                       /* ���������� �������� */
    double p_ratio, fg, q;          /* ��������������� ���������� */

    r = v_ncons[R];
    v = v_ncons[V];
    p = v_ncons[P];
    g = params->g;
    
    p_ratio = curr_press / p;
    if ( curr_press <= p ) {
        /* ����� ���������� */
        fg = 2.0 / ( g - 1.0 );
        *F = fg * c * ( pow( p_ratio, 1.0 / fg / g ) - 1.0 );
        *DF = ( 1.0 / r / c ) * pow( p_ratio, - 0.5 * ( g + 1.0 ) / g );
    }
    else {
        /* ������� ����� */
        q = sqrt( 0.5 * ( g + 1.0 ) / g * p_ratio + 0.5 * ( g - 1.0 ) / g );
        *F = ( curr_press - p ) / c / r / q;
        *DF = 0.25 * ( ( g + 1.0 ) * p_ratio + 3 * g - 1.0 ) / g / r / c / pow( q, 3.0 );
    }

}

/* ����������� ���������� ����������� ��� ������� �������� �� ���������� �������

   ���: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 157. - Subroutine GUESSP.

   params - ��������� � ����������� ��������������� ������������ (in)
   v_ncons_l[M] - ������ ����������� ���������� ����� �� ������� (in)
   v_ncons_r[M] - ������ ����������� ���������� ������ �� ������� (in)
   cl - �������� ����� ����� �� ������� (in)
   cr - �������� ����� ������ �� ������� (in)

   ���������� ������� ��������� ����������� */
double pressure_initial_guess( struct Parameters *params, double v_ncons_l[M], double v_ncons_r[M], double cl, double cr ) {

    double rl, vl, pl;                  /* ����������� ���������� ����� �� ������� */
    double rr, vr, pr;                  /* ����������� ���������� ������ �� ������� */
    double g;                           /* ���������� �������� */
    /* ��������� �����������, ������������ �� ����������� ������������ ��������������� �������
       � ����������� ���������� */
    double p_lin;
    double p_min, p_max;                /* ����������� � ������������ �������� ����� � ������ �� ������� */
    double p_ratio;                     /* ������� �� �������� ����� � ������ �� ������� */
    double p1, p2, g1, g2;              /* ��������������� ���������� ��� ������������� �������� */
    
    rl = v_ncons_l[R];
    vl = v_ncons_l[V];
    pl = v_ncons_l[P];
    rr = v_ncons_r[R];
    vr = v_ncons_r[V];
    pr = v_ncons_r[P];
    g = params->g;

    /* ��������� ����������� �� �������� ������
       Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
       1999. - P. 128. - Formula (4.47). */
    p_lin = max( 0.0, 0.5 * ( pl + pr ) - 0.125 * ( vr - vl ) * ( rl + rr ) * ( cl + cr ) );
    p_min = min( pl, pr );
    p_max = max( pl, pr );
    p_ratio = p_max / p_min;

    if ( ( p_ratio <= P_MAX_RATIO ) &&
       ( ( p_min < p_lin && p_lin < p_max ) || ( fabs( p_min - p_lin ) < EPS || fabs( p_max - p_lin ) < EPS ) ) ) {
        /* ��������� ����������� �� ��������������� ������ */
        return p_lin;
    } else {
        if ( p_lin < p_min ) {
            /* ��������� ����������� �� ���� ������ ����������
               Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
               1999. - P. 302. - Formula (9.32) + �������� �� ���������� ��������� ��������� */
            g1 = 0.5 * ( g - 1.0 ) / g;
            return pow( ( ( cl + cr - 0.5 * ( g - 1.0 ) * ( vr - vl ) ) / ( cl / pow( pl, g1 ) + cr / pow( pr, g1 ) ) ), 1.0 / g1 );
        } else {
            /* ��������� ����������� �� ���� ������� ������
               Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
               1999. - P. 128. - Formula (4.48) + �������� �� ���������� ��������� ��������� */
            g1 = 2.0 / ( g + 1.0 );
            g2 = ( g - 1.0 ) / ( g + 1.0 );
            p1 = sqrt( g1 / rl / ( g2 * pl + p_lin ) );
            p2 = sqrt( g1 / rr / ( g2 * pr + p_lin ) );
            return ( p1 * pl + p2 * pr - ( vr - vl ) ) / ( p1 + p2 );
        }
    }
    
}

/* ������� ������ �������

   ���: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 158. - Subroutine SAMPLE.

   params - ��������� � ����������� ��������������� ������������ (in)
   v_ncons_l[M] - ������ ���������������� ���������� ����� �� ������� (in)
   v_ncons_r[M] - ������ ���������������� ���������� ������ �� ������� (in)
   cl - �������� ����� ����� �� ������� (in)
   cr - �������� ����� ������ �� ������� (in)
   p_cont - �������� �� ���������� ������� (in)
   v_cont - �������� �� ���������� ������� (in)
   s - �������� x/t, ��� �������� ���������� ������� (in)

   v_ncons_res[M] - ������ ���������������� ���������� � ����������� ������������ (out) */
void sample_solid_solution( struct Parameters *params, double v_ncons_l[M], double v_ncons_r[M], double cl, double cr,
                            double p_cont, double v_cont, double s, double v_ncons_res[M] ) {

    double rl, vl, pl;                      /* ����������� ���������� ����� �� ������� */
    double rr, vr, pr;                      /* ����������� ���������� ������ �� ������� */
    double g1, g2, g3, g4, g5, g6, g7;      /* ��������������� ����������, ����������� �� ���������� ��������,
                                               � ������������ � Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer, 1999. - P. 153. */

    /* �������� ����� ���� */
    double shl, stl;        /* �������� "������" � "������" ����� ����� ���������� */
    double sl;              /* �������� ����� ������� ����� */

    /* �������� ������ ���� */
    double shr, str;        /* �������� "������" � "������" ������ ����� ���������� */
    double sr;              /* �������� ������ ������� ����� */

    double cml, cmr;        /* �������� ����� ����� � ������ �� ����������� ������� */
    double c;               /* ��������� �������� ����� ������ ����� ���������� */
    double p_ratio;
    double r, v, p;         /* ���������� �������� �������� ����, ���������, �������� � �������� */

    /* ��������������� ���������� */
    /* ��������� ����� �� ������� */
    rl = v_ncons_l[R];
    vl = v_ncons_l[V];
    pl = v_ncons_l[P];
    /* ��������� ������ �� ������� */
    rr = v_ncons_r[R];
    vr = v_ncons_r[V];
    pr = v_ncons_r[P];
    /* ����������� �� ���������� �������� */
    g1 = 0.5 * ( params->g - 1.0 ) / params->g;
    g2 = 0.5 * ( params->g + 1.0 ) / params->g;
    g3 = 2.0 * params->g / ( params->g - 1.0 );
    g4 = 2.0 / ( params->g - 1.0 );
    g5 = 2.0 / ( params->g + 1.0 );
    g6 = ( params->g - 1.0 ) / ( params->g + 1.0 );
    g7 = 0.5 * ( params->g - 1.0 );

    if ( s <= v_cont ) {
        /* ��������������� ����� - ����� �� ����������� ������� */
        if ( p_cont <= pl ) {
            /* ����� ����� ���������� */
            shl = vl - cl;
            if ( s <= shl ) {
                /* ��������� ����� �� ������� */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                cml = cl * pow( p_cont / pl, g1 );
                stl = v_cont - cml;
                if ( s > stl ) {
                    /* ��������� ����� �� ����������� ������� */
                    r = rl * pow( p_cont / pl, 1.0 / params->g );
                    v = v_cont;
                    p = p_cont;
                }
                else {
                    /* ��������� ������ ����� ����� ���������� */
                    v = g5 * ( cl + g7 * vl + s );
                    c = g5 * ( cl + g7 * ( vl - s ) );
                    r = rl * pow( c / cl, g4 );
                    p = pl * pow( c / cl, g3 );
                }
            }
        }
        else {
            /* ����� ������� ����� */
            p_ratio = p_cont / pl;
            sl = vl - cl * sqrt( g2 * p_ratio + g1 );
            if ( s <= sl ) {
                /* ��������� ����� �� ������� */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                /* ��������� �� ����� ������� ������ */
                r = rl * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
                v = v_cont;
                p = p_cont;
            }
        }
    }
    else {
        /* ��������������� ����� - ������ �� ����������� ������� */
        if ( p_cont > pr ) {
            /* ������ ������� ����� */
            p_ratio = p_cont / pr;
            sr = vr + cr * sqrt( g2 * p_ratio + g1 );
            if ( s >= sr ) {
                /* ��������� ������ �� ������� */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
                /* ��������� �� ������ ������� ������ */
                r = rr * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
                v = v_cont;
                p = p_cont;
            }
        }
        else {
            /* ������ ����� ���������� */
            shr = vr + cr;
            if ( s >= shr ) {
                /* ��������� ������ �� ������� */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
               cmr = cr * pow( p_cont / pr, g1 );
               str = v_cont + cmr;
               if ( s <= str ) {
                   /* ��������� ������ �� ����������� ������� */
                   r = rr * pow( p_cont / pr, 1.0 / params->g );
                   v = v_cont;
                   p = p_cont;
               }
               else {
                    /* ��������� ������ ������ ����� ���������� */
                    v = g5 * ( - cr + g7 * vr + s );
                    c = g5 * ( cr - g7 * ( vr - s ) );
                    r = rr * pow( c / cr, g4 );
                    p = pr * pow( c / cr, g3 );
               }
            }
        }
    }
    
    /* ������������ ��������� ������� � ����������� */
    v_ncons_res[R] = r;
    v_ncons_res[V] = v;
    v_ncons_res[P] = p;
    
}
