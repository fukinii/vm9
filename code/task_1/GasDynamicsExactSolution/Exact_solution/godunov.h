/*
 * godunov.h
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

#ifndef __GODUNOV_H_
#define __GODUNOV_H_

#include "main.h"

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
                                     double cl, double cr, double *p_cont, double *v_cont );

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
void calc_F_and_DF( struct Parameters *params, double curr_press, double v_ncons[M], double c, double *F, double *DF );

/* ����������� ���������� ����������� ��� ������� �������� �� ���������� �������

   ���: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 157. - Subroutine GUESSP.

   params - ��������� � ����������� ��������������� ������������ (in)
   v_ncons_l[M] - ������ ����������� ���������� ����� �� ������� (in)
   v_ncons_r[M] - ������ ����������� ���������� ������ �� ������� (in)
   cl - �������� ����� ����� �� ������� (in)
   cr - �������� ����� ������ �� ������� (in)

   ���������� ������� ��������� ����������� */
double pressure_initial_guess( struct Parameters *params, double v_ncons_l[M], double v_ncons_r[M], double cl, double cr );

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
                            double p_cont, double v_cont, double s, double v_ncons_res[M] );

#endif /* __GODUNOV_H_ */