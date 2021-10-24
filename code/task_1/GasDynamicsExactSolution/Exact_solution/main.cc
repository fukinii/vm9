/*
 * main.cc
 *
 * Построение точного решения задачи о распаде произвольного разрыва.
 *
 * (c) Уткин Павел, 2013
 *
 * Создан: 24 декабря 2013 г.
 *
 */

#include "main.h"

#include "io.h"
#include "grid.h"
#include "memory.h"
#include "utils.h"
#include "godunov.h"

int main( void ) {

    struct Parameters params;   /* структура с параметрами вычислительного эксперимента  */
    double *xc;                 /* массив координат центров ячеек сетки */
    double *x;                  /* массив координат узлов сетки */
    double cl, cr;              /* скорости звука слева и справа от разрыва */
    double p_cont, v_cont;      /* давление и скорость на контактном разрыве */
    int i_cell;                 /* индекс ячеек */
    int i_comp;                 /* индекс компонента вектора */
    double v_ncons_res[M];      /* текущее отобранное решение */
    double s;                   /* текущее значение автомодельной переменной */
    FILE *ex_sol_out;           /* файловый дескриптор для записи точного решения */
       
    printf( "\nExact solution of the Riemann problem for the gas dynamics equations\n(c) Pavel Utkin, ICAD RAS, MIPT, 2013\ne-mail: pavel_utk@mail.ru\n" );

    /* открытие файла для записи точного решения */
    fopen_s( &ex_sol_out, "exact_solution.dat", "wt" );
    if ( NULL == ex_sol_out ) {
        printf( "\nCan't open file exact_solution.dat for writing\n" );
    }

    /* считывание файла с параметрами задачи */
    read_parameters( &params );
 
    /* выделение памяти под массивы */
    get_memory_for_1D_double_array( params.cells_number, &xc );
    get_memory_for_1D_double_array( params.cells_number + 1, &x );

    /* определение координат центров ячеек сетки */
    build_grid( -0.5, 0.5, params.cells_number, xc, x );

    /* расчет скоростей звука слева и справа от разрыва */
    calc_sound_velocity( &params, params.left_params, &cl );
    calc_sound_velocity( &params, params.right_params, &cr );

    /* расчет давления и скорости на контактном разрыве */
    calc_contact_pressure_velocity( &params, params.left_params, params.right_params, cl, cr, &p_cont, &v_cont );

    /* цикл по ячейкам */
    for ( i_cell = 0; i_cell < params.cells_number; i_cell++ ) {
        /* изначально решение строится на отрезке [-0.5;0.5] */
        s = xc[i_cell] / params.stop_time;
        /* отбор решения для заданного значения s */
        sample_solid_solution( &params, params.left_params, params.right_params, cl, cr, p_cont, v_cont, s, v_ncons_res );
        /* запись вектора в файл */
        fprintf( ex_sol_out, "%e ", xc[i_cell] + 0.5 );
        for ( i_comp = 0; i_comp < M; i_comp++ )
            fprintf( ex_sol_out, "%e ", v_ncons_res[i_comp] );
        fprintf( ex_sol_out, "\n" );
    }

    /* освобождение памяти */
    free( xc );
    free( x );

    fclose( ex_sol_out );

    printf( "\nNormal finish. The exact solution is in exact_solution.dat in the current directory.\n\n" );

    return 0;

}