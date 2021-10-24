/*
 * utils.h
 *
 * Рабочие функции
 *
 * (c) Уткин Павел, 2012
 *
 * Создан: 24 марта 2012 г.
 *
 */

#ifndef __UTILS_H_
#define __UTILS_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "main.h"

/* Расчет скорости звука

   params - структура с параметрами вычислительного эксперимента (in)
   v_ncons[M] - вектор примитивных переменных (in)

   c - скорость звука (out) */
void calc_sound_velocity( struct Parameters *params, double v_ncons[M], double *c );

/* Функция max, возвращающая наибольшее из двух чисел
   
   a - первое число (in)
   b - второе число (in) */
double max( double a, double b );

/* Функция min, возвращающая наименьшее из двух чисел
   
   a - первое число (in)
   b - второе число (in) */
double min( double a, double b );

#endif /* __UTILS_H_ */