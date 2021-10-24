/*
 * main.cc
 *
 * ���������� ������� ������� ������ � ������� ������������� �������.
 *
 * (c) ����� �����, 2013
 *
 * ������: 24 ������� 2013 �.
 *
 */

#include "main.h"

#include "io.h"
#include "grid.h"
#include "memory.h"
#include "utils.h"
#include "godunov.h"

int main( void ) {

    struct Parameters params;   /* ��������� � ����������� ��������������� ������������  */
    double *xc;                 /* ������ ��������� ������� ����� ����� */
    double *x;                  /* ������ ��������� ����� ����� */
    double cl, cr;              /* �������� ����� ����� � ������ �� ������� */
    double p_cont, v_cont;      /* �������� � �������� �� ���������� ������� */
    int i_cell;                 /* ������ ����� */
    int i_comp;                 /* ������ ���������� ������� */
    double v_ncons_res[M];      /* ������� ���������� ������� */
    double s;                   /* ������� �������� ������������� ���������� */
    FILE *ex_sol_out;           /* �������� ���������� ��� ������ ������� ������� */
       
    printf( "\nExact solution of the Riemann problem for the gas dynamics equations\n(c) Pavel Utkin, ICAD RAS, MIPT, 2013\ne-mail: pavel_utk@mail.ru\n" );

    /* �������� ����� ��� ������ ������� ������� */
    fopen_s( &ex_sol_out, "exact_solution.dat", "wt" );
    if ( NULL == ex_sol_out ) {
        printf( "\nCan't open file exact_solution.dat for writing\n" );
    }

    /* ���������� ����� � ����������� ������ */
    read_parameters( &params );
 
    /* ��������� ������ ��� ������� */
    get_memory_for_1D_double_array( params.cells_number, &xc );
    get_memory_for_1D_double_array( params.cells_number + 1, &x );

    /* ����������� ��������� ������� ����� ����� */
    build_grid( -0.5, 0.5, params.cells_number, xc, x );

    /* ������ ��������� ����� ����� � ������ �� ������� */
    calc_sound_velocity( &params, params.left_params, &cl );
    calc_sound_velocity( &params, params.right_params, &cr );

    /* ������ �������� � �������� �� ���������� ������� */
    calc_contact_pressure_velocity( &params, params.left_params, params.right_params, cl, cr, &p_cont, &v_cont );

    /* ���� �� ������� */
    for ( i_cell = 0; i_cell < params.cells_number; i_cell++ ) {
        /* ���������� ������� �������� �� ������� [-0.5;0.5] */
        s = xc[i_cell] / params.stop_time;
        /* ����� ������� ��� ��������� �������� s */
        sample_solid_solution( &params, params.left_params, params.right_params, cl, cr, p_cont, v_cont, s, v_ncons_res );
        /* ������ ������� � ���� */
        fprintf( ex_sol_out, "%e ", xc[i_cell] + 0.5 );
        for ( i_comp = 0; i_comp < M; i_comp++ )
            fprintf( ex_sol_out, "%e ", v_ncons_res[i_comp] );
        fprintf( ex_sol_out, "\n" );
    }

    /* ������������ ������ */
    free( xc );
    free( x );

    fclose( ex_sol_out );

    printf( "\nNormal finish. The exact solution is in exact_solution.dat in the current directory.\n\n" );

    return 0;

}