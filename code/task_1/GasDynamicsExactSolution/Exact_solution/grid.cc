/*
 * grid.cc
 *
 * ������� ��� ������ � ��������� ������.
 *
 * (c) ����� �����, 2013
 *
 * ������: 17 ��� 2013 �.
 *
 */

#include "grid.h"

/* ����������� ��������� ����� � ������� ����� �����

   left_boundary_x - ���������� ����� ������� ��������� ������� (in)
   right_boundary_x - ���������� ������ ������� ��������� ������� (in)
   cells_num - ���������� ����� ����� (in)

   *xc - ������ ��������� ������� ����� (out)
   *x - ������ ��������� ����� ����� (out) */
void build_grid( double left_boundary_x, double right_boundary_x, int cells_num, double *xc, double *x ) {

    int i;
    double h;   /* ��� ����� */

   /* ���� ����������� ������ ����������� ������������� ����� ����� */
   h = ( right_boundary_x - left_boundary_x ) / cells_num;
   
   /* ���������� ����� */
   for ( i = 0; i < cells_num + 1; i++ ) {
       x[i] = left_boundary_x + i * h;
   }

   /* ���������� ������� ����� */
   for ( i = 0; i < cells_num; i++ ) {
       xc[i] = 0.5 * ( x[i] + x[i+1] );
   }

}
