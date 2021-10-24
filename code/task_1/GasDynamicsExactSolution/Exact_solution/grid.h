/*
 * grid.h
 *
 * ������� ��� ������ � ��������� ������.
 *
 * (c) ����� �����, 2013
 *
 * ������: 17 ��� 2013 �.
 *
 */

#ifndef __GRID_H_
#define __GRID_H_

#include "main.h"

/* ����������� ��������� ����� � ������� ����� �����

   left_boundary_x - ���������� ����� ������� ��������� ������� (in)
   right_boundary_x - ���������� ������ ������� ��������� ������� (in)
   cells_num - ���������� ����� ����� (in)

   *xc - ������ ��������� ������� ����� (out)
   *x - ������ ��������� ����� ����� (out) */
void build_grid( double left_boundary_x, double right_boundary_x, int cells_num, double *xc, double *x );

#endif /* __GRID_H_ */