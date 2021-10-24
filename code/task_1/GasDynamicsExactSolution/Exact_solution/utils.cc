/*
 * utils.cc
 *
 * Рабочие функции
 *
 * (c) Уткин Павел, 2012
 *
 * Создан: 24 марта 2012 г.
 *
 */

#include "utils.h"

/* Расчет скорости звука

   params - структура с параметрами вычислительного эксперимента (in)
   v_ncons[M] - вектор примитивных переменных (in)

   c - скорость звука (out) */
void calc_sound_velocity( struct Parameters *params, double v_ncons[M], double *c ) {

    *c = sqrt( params->g * v_ncons[P] / v_ncons[R] );

}

/* Функция max, возвращающая наибольшее из двух чисел
   
   a - первое число (in)
   b - второе число (in) */
double max( double a, double b ) {

    return ( a < b ) ? b : a;

}

/* Функция min, возвращающая наименьшее из двух чисел
   
   a - первое число (in)
   b - второе число (in) */
double min( double a, double b ) {

    return ( a < b ) ? a : b;

}