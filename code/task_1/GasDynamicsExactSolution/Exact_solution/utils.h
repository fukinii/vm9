/*
 * utils.h
 *
 * ������� �������
 *
 * (c) ����� �����, 2012
 *
 * ������: 24 ����� 2012 �.
 *
 */

#ifndef __UTILS_H_
#define __UTILS_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "main.h"

/* ������ �������� �����

   params - ��������� � ����������� ��������������� ������������ (in)
   v_ncons[M] - ������ ����������� ���������� (in)

   c - �������� ����� (out) */
void calc_sound_velocity( struct Parameters *params, double v_ncons[M], double *c );

/* ������� max, ������������ ���������� �� ���� �����
   
   a - ������ ����� (in)
   b - ������ ����� (in) */
double max( double a, double b );

/* ������� min, ������������ ���������� �� ���� �����
   
   a - ������ ����� (in)
   b - ������ ����� (in) */
double min( double a, double b );

#endif /* __UTILS_H_ */