/*
 * memory.cc
 *
 * ������� ��� ������ � �������.
 *
 * (c) ����� �����, 2013
 *
 * ������: 12 ������� 2013 �.
 *
 */

#include "memory.h"

/* ��������� ������ ��� ���������� ������ ��������� ���� double

   size - ������ ������� (in) 

   **a - ��������� �� ������ (out) */
void get_memory_for_1D_double_array( int size, double **a ) {

    *a = ( double * )malloc( size * sizeof( double ) );
    if ( NULL == *a ) {
        printf( "get_memory_for_1D_double_array -> Can't allocate memory.\n" );
        exit( EXIT_FAILURE );
    }

}

/* ��������� ������ ��� ��������� ������ ��������� ���� double

   size1 - ������ ����������� ������� (in)
   size2 - ������ ����������� ������� (in)

   ***a - ��������� �� ������ (out) */
void get_memory_for_2D_double_array( int size1, int size2, double ***a ) {

    int i;

    *a = ( double ** )malloc( size1 * sizeof( double * ) );
    if ( NULL == *a ) {
        printf( "get_memory_for_2D_double_array -> Can't allocate memory.\n" );
        exit( EXIT_FAILURE );
    }
    for ( i = 0; i < size1; i++ )
        get_memory_for_1D_double_array( size2, &((*a)[i]) );

}