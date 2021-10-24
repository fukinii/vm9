/*
 * io.cc
 *
 * Функции ввода/вывода.
 *
 * (c) Уткин Павел, 2013
 *
 * Создан: 17 мая 2012 г.
 *
 */

#include <stdio.h>

#include "io.h"

/* Считывание файла с параметрами задачи, заполнение полей соответствующей структуры и проверка
   корректности задания параметров

   *params - структура с параметрами вычислительного эксперимента (out) */
void read_parameters( struct Parameters *params ) {

    FILE *parameters;
    char string[MAX_STRING_SIZE];   /* для считывания строковой информации из файла */
        
    /* все параметры задачи находятся в файле parameters.txt */
    strcpy_s( string, "parameters.dat" );
    if ( ( fopen_s( &parameters, string, "rt" ) ) != 0 ) {
        printf( "\nread_parameters -> Can't open file %s for reading.\n\n", string );
        exit( EXIT_FAILURE );
    }

    /* считывание заголовка файла */
    fscanf_s( parameters, "%s", string, MAX_STRING_SIZE );

    /* сетка для построения точного решения */
    fscanf_s( parameters, "%s", string, MAX_STRING_SIZE );                                  /* заголовок раздела */
    fscanf_s( parameters, "%s %d", string, MAX_STRING_SIZE, &(params->cells_number) );      /* количество ячеек */
    if ( params->cells_number < 0 ) {
        printf( "\nread_parameters -> cells_number should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
   
    /* момент времени, для которого строится точное решение */
    fscanf_s( parameters, "%s", string, MAX_STRING_SIZE );                              /* заголовок раздела */
    fscanf_s( parameters, "%s %lf", string, MAX_STRING_SIZE, &(params->stop_time) );    /* момент времени, на который выводится результат */
    if ( params->stop_time <= 0.0 ) {
        printf( "\nread_parameters -> stop_time should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }

    /* начальные условия */
    fscanf_s( parameters, "%s", string, MAX_STRING_SIZE );                                  /* заголовок раздела */
    fscanf_s( parameters, "%s", string, MAX_STRING_SIZE );                                  /* заголовок - параметры слева от разрыва */
    fscanf_s( parameters, "%s %lf", string, MAX_STRING_SIZE, &(params->left_params[R]) );
    if ( params->left_params[R] <= 0.0 ) {
        printf( "\nread_parameters -> params->left_params[%d] should be positive.\n\n", R );
        exit( EXIT_FAILURE );
    }
    fscanf_s( parameters, "%s %lf", string, MAX_STRING_SIZE, &(params->left_params[V]) );
    fscanf_s( parameters, "%s %lf", string, MAX_STRING_SIZE, &(params->left_params[P]) );
    if ( params->left_params[P] <= 0.0 ) {
        printf( "\nread_parameters -> params->left_params[%d] should be positive.\n\n", P );
        exit( EXIT_FAILURE );
    }
    fscanf_s( parameters, "%s", string, MAX_STRING_SIZE );                                  /* заголовок - параметры справа от разрыва */
    fscanf_s( parameters, "%s %lf", string, MAX_STRING_SIZE, &(params->right_params[R]) );
    if ( params->right_params[R] <= 0.0 ) {
        printf( "\nread_parameters -> params->right_params[%d] should be positive.\n\n", R );
        exit( EXIT_FAILURE );
    }
    fscanf_s( parameters, "%s %lf", string, MAX_STRING_SIZE, &(params->right_params[V]) );
    fscanf_s( parameters, "%s %lf", string, MAX_STRING_SIZE, &(params->right_params[P]) );
    if ( params->right_params[P] <= 0.0 ) {
        printf( "\nread_parameters -> params->right_params[%d] should be positive.\n\n", P );
        exit( EXIT_FAILURE );
    }
        
    /* показатель адиабаты */
    fscanf_s( parameters, "%s", string, MAX_STRING_SIZE );                      /* заголовок раздела */
    fscanf_s( parameters, "%s %lf", string, MAX_STRING_SIZE, &(params->g) );    /* показатель адиабаты */
    if ( params->g <= 1.0 ) {
        printf( "\nread_parameters -> params->g should be grater than 1.0.\n\n" );
        exit( EXIT_FAILURE );
    }
    
}