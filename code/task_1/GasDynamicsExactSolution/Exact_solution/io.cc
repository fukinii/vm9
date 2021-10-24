/*
 * io.cc
 *
 * ������� �����/������.
 *
 * (c) ����� �����, 2013
 *
 * ������: 17 ��� 2012 �.
 *
 */

#include <stdio.h>

#include "io.h"

/* ���������� ����� � ����������� ������, ���������� ����� ��������������� ��������� � ��������
   ������������ ������� ����������

   *params - ��������� � ����������� ��������������� ������������ (out) */
void read_parameters( struct Parameters *params ) {

    FILE *parameters;
    char string[MAX_STRING_SIZE];   /* ��� ���������� ��������� ���������� �� ����� */
        
    /* ��� ��������� ������ ��������� � ����� parameters.txt */
    strcpy_s( string, "parameters.dat" );
    if ( ( fopen_s( &parameters, string, "rt" ) ) != 0 ) {
        printf( "\nread_parameters -> Can't open file %s for reading.\n\n", string );
        exit( EXIT_FAILURE );
    }

    /* ���������� ��������� ����� */
    fscanf_s( parameters, "%s", string, MAX_STRING_SIZE );

    /* ����� ��� ���������� ������� ������� */
    fscanf_s( parameters, "%s", string, MAX_STRING_SIZE );                                  /* ��������� ������� */
    fscanf_s( parameters, "%s %d", string, MAX_STRING_SIZE, &(params->cells_number) );      /* ���������� ����� */
    if ( params->cells_number < 0 ) {
        printf( "\nread_parameters -> cells_number should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }
   
    /* ������ �������, ��� �������� �������� ������ ������� */
    fscanf_s( parameters, "%s", string, MAX_STRING_SIZE );                              /* ��������� ������� */
    fscanf_s( parameters, "%s %lf", string, MAX_STRING_SIZE, &(params->stop_time) );    /* ������ �������, �� ������� ��������� ��������� */
    if ( params->stop_time <= 0.0 ) {
        printf( "\nread_parameters -> stop_time should be a positive value.\n\n" );
        exit( EXIT_FAILURE );
    }

    /* ��������� ������� */
    fscanf_s( parameters, "%s", string, MAX_STRING_SIZE );                                  /* ��������� ������� */
    fscanf_s( parameters, "%s", string, MAX_STRING_SIZE );                                  /* ��������� - ��������� ����� �� ������� */
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
    fscanf_s( parameters, "%s", string, MAX_STRING_SIZE );                                  /* ��������� - ��������� ������ �� ������� */
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
        
    /* ���������� �������� */
    fscanf_s( parameters, "%s", string, MAX_STRING_SIZE );                      /* ��������� ������� */
    fscanf_s( parameters, "%s %lf", string, MAX_STRING_SIZE, &(params->g) );    /* ���������� �������� */
    if ( params->g <= 1.0 ) {
        printf( "\nread_parameters -> params->g should be grater than 1.0.\n\n" );
        exit( EXIT_FAILURE );
    }
    
}