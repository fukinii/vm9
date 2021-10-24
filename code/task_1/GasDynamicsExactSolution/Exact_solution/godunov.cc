/*
 * godunov.cc
 *
 * Решение задачи Римана для системы уравнений газовой динамики.
 *
 * Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 3rd Edition. - Springer,
 * 2009. - P. 152 - 162.
 * 
 * (c) Уткин Павел, 2013
 *
 * Создан: 26 декабря 2013 г.
 *
 */

#include "godunov.h"

#include "utils.h"

/* Итерационная процедура расчета давления и скорости на контактном разрыве

   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 3rd Edition. - Springer,
   2009. - P. 155. - Subroutine STARPU.

   params - структура с параметрами вычислительного эксперимента (in)
   v_ncons_l[M] - вектор примитивных переменных слева от разрыва (in)
   v_ncons_r[M] - вектор примитивных переменных справа от разрыва (in)
   cl - скорость звука слева от разрыва (in)
   cr - скорость звука справа от разрыва (in)

   p_cont - давление на контактном разрыве (out)
   v_cont - скорость на контактном разрыве (out) */
void calc_contact_pressure_velocity( struct Parameters *params, double v_ncons_l[M], double v_ncons_r[M],
                                     double cl, double cr, double *p_cont, double *v_cont ) {

    double vl, vr;      /* скорости слева и справа от разрыва */
    double p_old;       /* значение давления на предыдущей итерации */
    double fl, fr;      /* значения функций */
    double fld, frd;    /* значения производных */
    int iter_num = 0;   /* количество проведенных итераций */
    double criteria;    /* переменная для определения сходимости */
    double g;           /* показатель адиабаты */

    /* введение обозначений для удобства записи формул */
    vl = v_ncons_l[V];
    vr = v_ncons_r[V];
    g = params->g;

    if ( 2.0 * ( cl + cr ) / ( g - 1.0 ) <= vr - vl ) {
        /* случай возникновения вакуума */
        printf( "\ncalc_contact_pressure_velocity -> vacuum is generated\n" );
        exit( EXIT_FAILURE );
    }

    /* расчет начального приближения для давления */
    p_old = pressure_initial_guess( params, v_ncons_l, v_ncons_r, cl, cr );
    if ( p_old < 0.0 ) {
        printf( "\ncalc_contact_pressure_velocity -> initial pressure guess is negative " );
        exit( EXIT_FAILURE );
    }
    
    /* решение нелинейного уравнения для нахождения давления на контактном разрыве методом Ньютона-Рафсона */
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

    /* скорость контактного разрыва */
    *v_cont = 0.5 * ( vl + vr + fr - fl );

}

/* Расчет функции F, определяющей скорость газа на контактном разрыве, и ее производной по давлению среды DF

   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 3rd Edition. - Springer,
   1999. - P. 158. - Subroutine PREFUN.
   + поправка на двучленное уравнение состояния:

   params - структура с параметрами вычислительного эксперимента (in)
   curr_press - давление с предыдущей итерации (in)
   v_ncons[M] - вектор примитивных переменных (in)
   c - скорость звука в невозмущенной среде (in)

   F - значение функции (out)
   DF - значение производной (out) */
void calc_F_and_DF( struct Parameters *params, double curr_press, double v_ncons[M], double c, double *F, double *DF ) {

    double r, v, p;                 /* примитивные переменные */
    double g;                       /* показатель адиабаты */
    double p_ratio, fg, q;          /* вспомогательные переменные */

    r = v_ncons[R];
    v = v_ncons[V];
    p = v_ncons[P];
    g = params->g;
    
    p_ratio = curr_press / p;
    if ( curr_press <= p ) {
        /* волна разрежения */
        fg = 2.0 / ( g - 1.0 );
        *F = fg * c * ( pow( p_ratio, 1.0 / fg / g ) - 1.0 );
        *DF = ( 1.0 / r / c ) * pow( p_ratio, - 0.5 * ( g + 1.0 ) / g );
    }
    else {
        /* ударная волна */
        q = sqrt( 0.5 * ( g + 1.0 ) / g * p_ratio + 0.5 * ( g - 1.0 ) / g );
        *F = ( curr_press - p ) / c / r / q;
        *DF = 0.25 * ( ( g + 1.0 ) * p_ratio + 3 * g - 1.0 ) / g / r / c / pow( q, 3.0 );
    }

}

/* Определение начального приближения для расчета давления на контактном разрыве

   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 157. - Subroutine GUESSP.

   params - структура с параметрами вычислительного эксперимента (in)
   v_ncons_l[M] - вектор примитивных переменных слева от разрыва (in)
   v_ncons_r[M] - вектор примитивных переменных справа от разрыва (in)
   cl - скорость звука слева от разрыва (in)
   cr - скорость звука справа от разрыва (in)

   Возвращает искомое начальное приближения */
double pressure_initial_guess( struct Parameters *params, double v_ncons_l[M], double v_ncons_r[M], double cl, double cr ) {

    double rl, vl, pl;                  /* примитивные переменные слева от разрыва */
    double rr, vr, pr;                  /* примитивные переменные справа от разрыва */
    double g;                           /* показатель адиабаты */
    /* начальное приближение, рассчитанное на освановании рассмотрения линеаризованной системы
       в примитивных переменных */
    double p_lin;
    double p_min, p_max;                /* минимальное и максимальное давления слева и справа от разрыва */
    double p_ratio;                     /* перепад по давлению слева и справа от разрыва */
    double p1, p2, g1, g2;              /* вспомогательные переменные для промежуточных расчетов */
    
    rl = v_ncons_l[R];
    vl = v_ncons_l[V];
    pl = v_ncons_l[P];
    rr = v_ncons_r[R];
    vr = v_ncons_r[V];
    pr = v_ncons_r[P];
    g = params->g;

    /* Начальное приближение из линейной задачи
       Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
       1999. - P. 128. - Formula (4.47). */
    p_lin = max( 0.0, 0.5 * ( pl + pr ) - 0.125 * ( vr - vl ) * ( rl + rr ) * ( cl + cr ) );
    p_min = min( pl, pr );
    p_max = max( pl, pr );
    p_ratio = p_max / p_min;

    if ( ( p_ratio <= P_MAX_RATIO ) &&
       ( ( p_min < p_lin && p_lin < p_max ) || ( fabs( p_min - p_lin ) < EPS || fabs( p_max - p_lin ) < EPS ) ) ) {
        /* Начальное приближение из линеаризованной задачи */
        return p_lin;
    } else {
        if ( p_lin < p_min ) {
            /* Начальное приближение по двум волнам разрежения
               Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
               1999. - P. 302. - Formula (9.32) + поправка на двучленное уравнение состояния */
            g1 = 0.5 * ( g - 1.0 ) / g;
            return pow( ( ( cl + cr - 0.5 * ( g - 1.0 ) * ( vr - vl ) ) / ( cl / pow( pl, g1 ) + cr / pow( pr, g1 ) ) ), 1.0 / g1 );
        } else {
            /* Начальное приближение по двум ударным волнам
               Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
               1999. - P. 128. - Formula (4.48) + поправка на двучленное уравнение состояния */
            g1 = 2.0 / ( g + 1.0 );
            g2 = ( g - 1.0 ) / ( g + 1.0 );
            p1 = sqrt( g1 / rl / ( g2 * pl + p_lin ) );
            p2 = sqrt( g1 / rr / ( g2 * pr + p_lin ) );
            return ( p1 * pl + p2 * pr - ( vr - vl ) ) / ( p1 + p2 );
        }
    }
    
}

/* Функция отбора решения

   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 158. - Subroutine SAMPLE.

   params - структура с параметрами вычислительного эксперимента (in)
   v_ncons_l[M] - вектор неконсервативных переменных слева от разрыва (in)
   v_ncons_r[M] - вектор неконсервативных переменных справа от разрыва (in)
   cl - скорость звука слева от разрыва (in)
   cr - скорость звука справа от разрыва (in)
   p_cont - давление на контактном разрыве (in)
   v_cont - скорость на контактном разрыве (in)
   s - значение x/t, для которого отбирается решение (in)

   v_ncons_res[M] - вектор неконсервативных переменных с отобранными компонентами (out) */
void sample_solid_solution( struct Parameters *params, double v_ncons_l[M], double v_ncons_r[M], double cl, double cr,
                            double p_cont, double v_cont, double s, double v_ncons_res[M] ) {

    double rl, vl, pl;                      /* примитивные переменные слева от разрыва */
    double rr, vr, pr;                      /* примитивные переменные справа от разрыва */
    double g1, g2, g3, g4, g5, g6, g7;      /* вспомогательные переменные, производные от показателя адиабаты,
                                               в соответствии с Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer, 1999. - P. 153. */

    /* скорости левых волн */
    double shl, stl;        /* скорости "головы" и "хвоста" левой волны разрежения */
    double sl;              /* скорость левой ударной волны */

    /* скорости правых волн */
    double shr, str;        /* скорости "головы" и "хвоста" правой волны разрежения */
    double sr;              /* скорость правой ударной волны */

    double cml, cmr;        /* скорости звука слева и справа от контактного разрыва */
    double c;               /* локальная скорость звука внутри волны разрежения */
    double p_ratio;
    double r, v, p;         /* отобранные значения объемной доли, плотности, скорости и давления */

    /* вспомогательные переменные */
    /* параметры слева от разрыва */
    rl = v_ncons_l[R];
    vl = v_ncons_l[V];
    pl = v_ncons_l[P];
    /* параметры справа от разрыва */
    rr = v_ncons_r[R];
    vr = v_ncons_r[V];
    pr = v_ncons_r[P];
    /* производные от показателя адиабаты */
    g1 = 0.5 * ( params->g - 1.0 ) / params->g;
    g2 = 0.5 * ( params->g + 1.0 ) / params->g;
    g3 = 2.0 * params->g / ( params->g - 1.0 );
    g4 = 2.0 / ( params->g - 1.0 );
    g5 = 2.0 / ( params->g + 1.0 );
    g6 = ( params->g - 1.0 ) / ( params->g + 1.0 );
    g7 = 0.5 * ( params->g - 1.0 );

    if ( s <= v_cont ) {
        /* рассматриваемая точка - слева от контактного разрыва */
        if ( p_cont <= pl ) {
            /* левая волна разрежения */
            shl = vl - cl;
            if ( s <= shl ) {
                /* параметры слева от разрыва */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                cml = cl * pow( p_cont / pl, g1 );
                stl = v_cont - cml;
                if ( s > stl ) {
                    /* параметры слева от контактного разрыва */
                    r = rl * pow( p_cont / pl, 1.0 / params->g );
                    v = v_cont;
                    p = p_cont;
                }
                else {
                    /* параметры внутри левой волны разрежения */
                    v = g5 * ( cl + g7 * vl + s );
                    c = g5 * ( cl + g7 * ( vl - s ) );
                    r = rl * pow( c / cl, g4 );
                    p = pl * pow( c / cl, g3 );
                }
            }
        }
        else {
            /* левая ударная волна */
            p_ratio = p_cont / pl;
            sl = vl - cl * sqrt( g2 * p_ratio + g1 );
            if ( s <= sl ) {
                /* параметры слева от разрыва */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                /* параметры за левой ударной волной */
                r = rl * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
                v = v_cont;
                p = p_cont;
            }
        }
    }
    else {
        /* рассматриваемая точка - справа от контактного разрыва */
        if ( p_cont > pr ) {
            /* правая ударная волна */
            p_ratio = p_cont / pr;
            sr = vr + cr * sqrt( g2 * p_ratio + g1 );
            if ( s >= sr ) {
                /* параметры справа от разрыва */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
                /* параметры за правой ударной волной */
                r = rr * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
                v = v_cont;
                p = p_cont;
            }
        }
        else {
            /* правая волна разрежения */
            shr = vr + cr;
            if ( s >= shr ) {
                /* параметры справа от разрыва */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
               cmr = cr * pow( p_cont / pr, g1 );
               str = v_cont + cmr;
               if ( s <= str ) {
                   /* параметры справа от контактного разрыва */
                   r = rr * pow( p_cont / pr, 1.0 / params->g );
                   v = v_cont;
                   p = p_cont;
               }
               else {
                    /* параметры внутри правой волны разрежения */
                    v = g5 * ( - cr + g7 * vr + s );
                    c = g5 * ( cr - g7 * ( vr - s ) );
                    r = rr * pow( c / cr, g4 );
                    p = pr * pow( c / cr, g3 );
               }
            }
        }
    }
    
    /* формирование выходного вектора с результатом */
    v_ncons_res[R] = r;
    v_ncons_res[V] = v;
    v_ncons_res[P] = p;
    
}
