/*
 * godunov.h
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

#ifndef __GODUNOV_H_
#define __GODUNOV_H_

#include "main.h"

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
                                     double cl, double cr, double *p_cont, double *v_cont );

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
void calc_F_and_DF( struct Parameters *params, double curr_press, double v_ncons[M], double c, double *F, double *DF );

/* Определение начального приближения для расчета давления на контактном разрыве

   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 157. - Subroutine GUESSP.

   params - структура с параметрами вычислительного эксперимента (in)
   v_ncons_l[M] - вектор примитивных переменных слева от разрыва (in)
   v_ncons_r[M] - вектор примитивных переменных справа от разрыва (in)
   cl - скорость звука слева от разрыва (in)
   cr - скорость звука справа от разрыва (in)

   Возвращает искомое начальное приближения */
double pressure_initial_guess( struct Parameters *params, double v_ncons_l[M], double v_ncons_r[M], double cl, double cr );

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
                            double p_cont, double v_cont, double s, double v_ncons_res[M] );

#endif /* __GODUNOV_H_ */