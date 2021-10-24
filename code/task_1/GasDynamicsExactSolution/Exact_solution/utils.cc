/*
 * utils.cc
 *
 * ������� �������
 *
 * (c) ����� �����, 2012
 *
 * ������: 24 ����� 2012 �.
 *
 */

#include "utils.h"

/* ������ �������� �����

   params - ��������� � ����������� ��������������� ������������ (in)
   v_ncons[M] - ������ ����������� ���������� (in)

   c - �������� ����� (out) */
void calc_sound_velocity( struct Parameters *params, double v_ncons[M], double *c ) {

    *c = sqrt( params->g * v_ncons[P] / v_ncons[R] );

}

/* ������� max, ������������ ���������� �� ���� �����
   
   a - ������ ����� (in)
   b - ������ ����� (in) */
double max( double a, double b ) {

    return ( a < b ) ? b : a;

}

/* ������� min, ������������ ���������� �� ���� �����
   
   a - ������ ����� (in)
   b - ������ ����� (in) */
double min( double a, double b ) {

    return ( a < b ) ? a : b;

}