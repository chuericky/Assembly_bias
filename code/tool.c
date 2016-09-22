// Updated: Nov 18, 2015

#include "tool.h"
#include <stdio.h>
#include <stdlib.h>

FILE *open_check(char *filename, char *mode) {
    // Open file
    FILE *fp = fopen(filename, mode);
    if (!fp) {
        fprintf(stderr, "Error opening file %s with mode %s\n", filename, mode);
        perror("fopen"); exit(EXIT_FAILURE);
    }
    return fp;
}

int count_line (int spot, FILE *countfile, int line)
{
    // A function to count number of lines in a file.
    do {
        spot = fgetc (countfile);
        if (spot == '\n')
            line++;
    } while (spot != EOF);
    
    // last line doesn't end with a new line!
    // but there has to be a line at least before the last line
    if (spot != '\n' && line != 0)
        line++;
    
    if (line == 0)
        line++;
    
    rewind (countfile);
    
    return line;
}

int jackarea (char field[]) {
    // A function to determine the jackknife region for the galaxies.
    int field_ind = 0;
    
    // W1
    if (strcmp("W1p4p3", field) == 0 || strcmp("W1p3p3", field) == 0 || strcmp("W1p2p3", field) == 0 || strcmp("W1p1p3", field) == 0)
        field_ind = 0;
    if (strcmp("W1m0p3", field) == 0 || strcmp("W1m1p3", field) == 0 || strcmp("W1m2p3", field) == 0 || strcmp("W1m3p3", field) == 0)
        field_ind = 1;
    if (strcmp("W1p4p2", field) == 0 || strcmp("W1p3p2", field) == 0 || strcmp("W1p2p2", field) == 0 || strcmp("W1p1p2", field) == 0)
        field_ind = 2;
    if (strcmp("W1m0p2", field) == 0 || strcmp("W1m1p2", field) == 0 || strcmp("W1m2p2", field) == 0 || strcmp("W1m3p2", field) == 0)
        field_ind = 3;
    if (strcmp("W1p4p1", field) == 0 || strcmp("W1p3p1", field) == 0 || strcmp("W1p2p1", field) == 0 || strcmp("W1p1p1", field) == 0)
        field_ind = 4;
    if (strcmp("W1m0p1", field) == 0 || strcmp("W1m1p1", field) == 0 || strcmp("W1m2p1", field) == 0 || strcmp("W1m3p1", field) == 0)
        field_ind = 5;
    if (strcmp("W1p4m0", field) == 0 || strcmp("W1p3m0", field) == 0 || strcmp("W1p2m0", field) == 0 || strcmp("W1p1m0", field) == 0)
        field_ind = 6;
    if (strcmp("W1m0m0", field) == 0 || strcmp("W1m1m0", field) == 0 || strcmp("W1m2m0", field) == 0 || strcmp("W1m3m0", field) == 0)
        field_ind = 7;
    if (strcmp("W1p4m1", field) == 0 || strcmp("W1p3m1", field) == 0 || strcmp("W1p2m1", field) == 0 || strcmp("W1p1m1", field) == 0)
        field_ind = 8;
    if (strcmp("W1m0m1", field) == 0 || strcmp("W1m1m1", field) == 0 || strcmp("W1m2m1", field) == 0 || strcmp("W1m3m1", field) == 0)
        field_ind = 9;
    if (strcmp("W1p4m2", field) == 0 || strcmp("W1p3m2", field) == 0 || strcmp("W1p2m2", field) == 0 || strcmp("W1p1m2", field) == 0)
        field_ind = 10;
    if (strcmp("W1m0m2", field) == 0 || strcmp("W1m1m2", field) == 0 || strcmp("W1m2m2", field) == 0 || strcmp("W1m3m2", field) == 0)
        field_ind = 11;
    if (strcmp("W1p4m3", field) == 0 || strcmp("W1p3m3", field) == 0 || strcmp("W1p2m3", field) == 0 || strcmp("W1p1m3", field) == 0)
        field_ind = 12;
    if (strcmp("W1m0m3", field) == 0 || strcmp("W1m1m3", field) == 0 || strcmp("W1m2m3", field) == 0 || strcmp("W1m3m3", field) == 0)
        field_ind = 13;
    if (strcmp("W1p4m4", field) == 0 || strcmp("W1p3m4", field) == 0 || strcmp("W1p2m4", field) == 0 || strcmp("W1p1m4", field) == 0)
        field_ind = 14;
    if (strcmp("W1m0m4", field) == 0 || strcmp("W1m1m4", field) == 0 || strcmp("W1m2m4", field) == 0 || strcmp("W1m3m4", field) == 0)
        field_ind = 15;
    if (strcmp("W1m4p3", field) == 0 || strcmp("W1m4p2", field) == 0 || strcmp("W1m4p1", field) == 0 || strcmp("W1m4m0", field) == 0)
        field_ind = 16;
    if (strcmp("W1m4m1", field) == 0 || strcmp("W1m4m2", field) == 0 || strcmp("W1m4m3", field) == 0 || strcmp("W1m4m4", field) == 0)
        field_ind = 17;
    
    // W2
    if (strcmp("W2p3p2", field) == 0 || strcmp("W2p3p1", field) == 0 || strcmp("W2p3m0", field) == 0 || strcmp("W2p3m1", field) == 0)
        field_ind = 18;
    if (strcmp("W2p2p2", field) == 0 || strcmp("W2p2p1", field) == 0 || strcmp("W2p2m0", field) == 0 || strcmp("W2p2m1", field) == 0)
        field_ind = 19;
    if (strcmp("W2p1p2", field) == 0 || strcmp("W2p1p1", field) == 0 || strcmp("W2p1m0", field) == 0 || strcmp("W2p1m1", field) == 0)
        field_ind = 20;
    if (strcmp("W2m0p2", field) == 0 || strcmp("W2m0p1", field) == 0 || strcmp("W2m0m0", field) == 0 || strcmp("W2m0m1", field) == 0)
        field_ind = 21;
    if (strcmp("W2m1p2", field) == 0 || strcmp("W2m1p1", field) == 0 || strcmp("W2m1m0", field) == 0 || strcmp("W2m1m1", field) == 0)
        field_ind = 22;
    if (strcmp("W2p3p3", field) == 0 || strcmp("W2p2p3", field) == 0 || strcmp("W2p1p3", field) == 0 || strcmp("W2m0p3", field) == 0 || strcmp("W2m1p3", field) == 0)
        field_ind = 23;
    
    // W3
    if (strcmp("W3p3p3", field) == 0 || strcmp("W3p2p3", field) == 0)
        field_ind = 24;
    if (strcmp("W3p1p3", field) == 0 || strcmp("W3m0p3", field) == 0)
        field_ind = 25;
    if (strcmp("W3m1p3", field) == 0 || strcmp("W3m2p3", field) == 0)
        field_ind = 26;
    if (strcmp("W3p3p2", field) == 0 || strcmp("W3p2p2", field) == 0)
        field_ind = 27;
    if (strcmp("W3p1p2", field) == 0 || strcmp("W3m0p2", field) == 0)
        field_ind = 28;
    if (strcmp("W3m1p2", field) == 0 || strcmp("W3m2p2", field) == 0)
        field_ind = 29;
    if (strcmp("W3p3p1", field) == 0 || strcmp("W3p2p1", field) == 0)
        field_ind = 30;
    if (strcmp("W3p1p1", field) == 0 || strcmp("W3m0p1", field) == 0)
        field_ind = 31;
    if (strcmp("W3m1p1", field) == 0 || strcmp("W3m2p1", field) == 0)
        field_ind = 32;
    if (strcmp("W3p3m0", field) == 0 || strcmp("W3p2m0", field) == 0)
        field_ind = 33;
    if (strcmp("W3p1m0", field) == 0 || strcmp("W3m0m0", field) == 0)
        field_ind = 34;
    if (strcmp("W3m1m0", field) == 0 || strcmp("W3m2m0", field) == 0)
        field_ind = 35;
    if (strcmp("W3p3m1", field) == 0 || strcmp("W3p2m1", field) == 0)
        field_ind = 36;
    if (strcmp("W3p1m1", field) == 0 || strcmp("W3m0m1", field) == 0)
        field_ind = 37;
    if (strcmp("W3m1m1", field) == 0 || strcmp("W3m2m1", field) == 0)
        field_ind = 38;
    if (strcmp("W3p3m2", field) == 0 || strcmp("W3p2m2", field) == 0)
        field_ind = 39;
    if (strcmp("W3p1m2", field) == 0 || strcmp("W3m0m2", field) == 0)
        field_ind = 40;
    if (strcmp("W3m1m2", field) == 0 || strcmp("W3m2m2", field) == 0)
        field_ind = 41;
    if (strcmp("W3p3m3", field) == 0 || strcmp("W3p2m3", field) == 0)
        field_ind = 42;
    if (strcmp("W3p1m3", field) == 0 || strcmp("W3m0m3", field) == 0)
        field_ind = 43;
    if (strcmp("W3m1m3", field) == 0 || strcmp("W3m2m3", field) == 0 || strcmp("W3m3m3", field) == 0)
        field_ind = 44;
    if (strcmp("W3m3m2", field) == 0 || strcmp("W3m3m1", field) == 0)
        field_ind = 45;
    if (strcmp("W3m3m0", field) == 0 || strcmp("W3m3p1", field) == 0)
        field_ind = 46;
    if (strcmp("W3m3p2", field) == 0 || strcmp("W3m3p3", field) == 0)
        field_ind = 47;
    
    // W4
    if (strcmp("W4p2m2", field) == 0 || strcmp("W4p1m2", field) == 0 || strcmp("W4m0m2", field) == 0 || strcmp("W4m1m2", field) == 0)
        field_ind = 48;
    if (strcmp("W4p2m1", field) == 0 || strcmp("W4p1m1", field) == 0 || strcmp("W4m0m1", field) == 0 || strcmp("W4m1m1", field) == 0)
        field_ind = 49;
    if (strcmp("W4p2m0", field) == 0 || strcmp("W4p1m0", field) == 0 || strcmp("W4m0m0", field) == 0 || strcmp("W4m1m0", field) == 0)
        field_ind = 50;
    if (strcmp("W4p1p1", field) == 0 || strcmp("W4m0p1", field) == 0 || strcmp("W4m1p1", field) == 0 || strcmp("W4m1p2", field) == 0 || strcmp("W4m1p3", field) == 0)
        field_ind = 51;
    if (strcmp("W4m2m0", field) == 0 || strcmp("W4m2p1", field) == 0 || strcmp("W4m2p2", field) == 0 || strcmp("W4m2p3", field) == 0)
        field_ind = 52;
    if (strcmp("W4m3m0", field) == 0 || strcmp("W4m3p1", field) == 0 || strcmp("W4m3p2", field) == 0 || strcmp("W4m3p3", field) == 0)
        field_ind = 53;
    
    return field_ind;
}

int jackarea2 (char field[]) {
    // A function to determine the jackknife region for the galaxies.
    int field_ind = 0;
    
    // W1
    if (strcmp("W1p4p3", field) == 0 || strcmp("W1p3p3", field) == 0)
        field_ind = 0;
    if (strcmp("W1p2p3", field) == 0 || strcmp("W1p1p3", field) == 0)
        field_ind = 1;
    if (strcmp("W1m0p3", field) == 0 || strcmp("W1m1p3", field) == 0)
        field_ind = 2;
    if (strcmp("W1m2p3", field) == 0 || strcmp("W1m3p3", field) == 0)
        field_ind = 3;
    if (strcmp("W1p4p2", field) == 0 || strcmp("W1p3p2", field) == 0)
        field_ind = 4;
    if (strcmp("W1p2p2", field) == 0 || strcmp("W1p1p2", field) == 0)
        field_ind = 5;
    if (strcmp("W1m0p2", field) == 0 || strcmp("W1m1p2", field) == 0)
        field_ind = 6;
    if (strcmp("W1m2p2", field) == 0 || strcmp("W1m3p2", field) == 0)
        field_ind = 7;
    if (strcmp("W1p4p1", field) == 0 || strcmp("W1p3p1", field) == 0)
        field_ind = 8;
    if (strcmp("W1p2p1", field) == 0 || strcmp("W1p1p1", field) == 0)
        field_ind = 9;
    if (strcmp("W1m0p1", field) == 0 || strcmp("W1m1p1", field) == 0)
        field_ind = 10;
    if (strcmp("W1m2p1", field) == 0 || strcmp("W1m3p1", field) == 0)
        field_ind = 11;
    if (strcmp("W1p4m0", field) == 0 || strcmp("W1p3m0", field) == 0)
        field_ind = 12;
    if (strcmp("W1p2m0", field) == 0 || strcmp("W1p1m0", field) == 0)
        field_ind = 13;
    if (strcmp("W1m0m0", field) == 0 || strcmp("W1m1m0", field) == 0)
        field_ind = 14;
    if (strcmp("W1m2m0", field) == 0 || strcmp("W1m3m0", field) == 0)
        field_ind = 15;
    if (strcmp("W1p4m1", field) == 0 || strcmp("W1p3m1", field) == 0)
        field_ind = 16;
    if (strcmp("W1p2m1", field) == 0 || strcmp("W1p1m1", field) == 0)
        field_ind = 17;
    if (strcmp("W1m0m1", field) == 0 || strcmp("W1m1m1", field) == 0)
        field_ind = 18;
    if (strcmp("W1m2m1", field) == 0 || strcmp("W1m3m1", field) == 0)
        field_ind = 19;
    if (strcmp("W1p4m2", field) == 0 || strcmp("W1p3m2", field) == 0)
        field_ind = 20;
    if (strcmp("W1p2m2", field) == 0 || strcmp("W1p1m2", field) == 0)
        field_ind = 21;
    if (strcmp("W1m0m2", field) == 0 || strcmp("W1m1m2", field) == 0)
        field_ind = 22;
    if (strcmp("W1m2m2", field) == 0 || strcmp("W1m3m2", field) == 0)
        field_ind = 23;
    if (strcmp("W1p4m3", field) == 0 || strcmp("W1p3m3", field) == 0)
        field_ind = 24;
    if (strcmp("W1p2m3", field) == 0 || strcmp("W1p1m3", field) == 0)
        field_ind = 25;
    if (strcmp("W1m0m3", field) == 0 || strcmp("W1m1m3", field) == 0)
        field_ind = 26;
    if (strcmp("W1m2m3", field) == 0 || strcmp("W1m3m3", field) == 0)
        field_ind = 27;
    if (strcmp("W1p4m4", field) == 0 || strcmp("W1p3m4", field) == 0)
        field_ind = 28;
    if (strcmp("W1p2m4", field) == 0 || strcmp("W1p1m4", field) == 0)
        field_ind = 29;
    if (strcmp("W1m0m4", field) == 0 || strcmp("W1m1m4", field) == 0)
        field_ind = 30;
    if (strcmp("W1m2m4", field) == 0 || strcmp("W1m3m4", field) == 0)
        field_ind = 31;
    if (strcmp("W1m4p3", field) == 0 || strcmp("W1m4p2", field) == 0)
        field_ind = 32;
    if (strcmp("W1m4p1", field) == 0 || strcmp("W1m4m0", field) == 0)
        field_ind = 33;
    if (strcmp("W1m4m1", field) == 0 || strcmp("W1m4m2", field) == 0)
        field_ind = 34;
    if (strcmp("W1m4m3", field) == 0 || strcmp("W1m4m4", field) == 0)
        field_ind = 35;
    
    // W2
    if (strcmp("W2p3p2", field) == 0 || strcmp("W2p3p1", field) == 0)
        field_ind = 36;
    if (strcmp("W2p3m0", field) == 0 || strcmp("W2p3m1", field) == 0)
        field_ind = 37;
    if (strcmp("W2p2p2", field) == 0 || strcmp("W2p2p1", field) == 0)
        field_ind = 38;
    if (strcmp("W2p2m0", field) == 0 || strcmp("W2p2m1", field) == 0)
        field_ind = 39;
    if (strcmp("W2p1p2", field) == 0 || strcmp("W2p1p1", field) == 0)
        field_ind = 40;
    if (strcmp("W2p1m0", field) == 0 || strcmp("W2p1m1", field) == 0)
        field_ind = 41;
    if (strcmp("W2m0p2", field) == 0 || strcmp("W2m0p1", field) == 0)
        field_ind = 42;
    if (strcmp("W2m0m0", field) == 0 || strcmp("W2m0m1", field) == 0)
        field_ind = 43;
    if (strcmp("W2m1p2", field) == 0 || strcmp("W2m1p1", field) == 0)
        field_ind = 44;
    if (strcmp("W2m1m0", field) == 0 || strcmp("W2m1m1", field) == 0)
        field_ind = 45;
    if (strcmp("W2p3p3", field) == 0 || strcmp("W2p2p3", field) == 0)
        field_ind = 46;
    if (strcmp("W2p1p3", field) == 0 || strcmp("W2m0p3", field) == 0 || strcmp("W2m1p3", field) == 0)
        field_ind = 47;
    
    // W3
    if (strcmp("W3p3p3", field) == 0)
        field_ind = 48;
    if (strcmp("W3p2p3", field) == 0)
        field_ind = 49;
    if (strcmp("W3p1p3", field) == 0)
    field_ind = 50;
    if (strcmp("W3m0p3", field) == 0)
    field_ind = 51;
    if (strcmp("W3m1p3", field) == 0)
    field_ind = 52;
    if (strcmp("W3m2p3", field) == 0)
    field_ind = 53;
    if (strcmp("W3p3p2", field) == 0)
    field_ind = 54;
    if (strcmp("W3p2p2", field) == 0)
    field_ind = 55;
    if (strcmp("W3p1p2", field) == 0)
    field_ind = 56;
    if (strcmp("W3m0p2", field) == 0)
    field_ind = 57;
    if (strcmp("W3m1p2", field) == 0)
    field_ind = 58;
    if (strcmp("W3m2p2", field) == 0)
    field_ind = 59;
    if (strcmp("W3p3p1", field) == 0)
    field_ind = 60;
    if (strcmp("W3p2p1", field) == 0)
    field_ind = 61;
    if (strcmp("W3p1p1", field) == 0)
    field_ind = 62;
    if (strcmp("W3m0p1", field) == 0)
    field_ind = 63;
    if (strcmp("W3m1p1", field) == 0)
    field_ind = 64;
    if (strcmp("W3m2p1", field) == 0)
    field_ind = 65;
    if (strcmp("W3p3m0", field) == 0)
    field_ind = 66;
    if (strcmp("W3p2m0", field) == 0)
    field_ind = 67;
    if (strcmp("W3p1m0", field) == 0)
    field_ind = 68;
    if (strcmp("W3m0m0", field) == 0)
    field_ind = 69;
    if (strcmp("W3m1m0", field) == 0)
    field_ind = 70;
    if (strcmp("W3m2m0", field) == 0)
    field_ind = 71;
    if (strcmp("W3p3m1", field) == 0)
    field_ind = 72;
    if (strcmp("W3p2m1", field) == 0)
    field_ind = 73;
    if (strcmp("W3p1m1", field) == 0)
    field_ind = 74;
    if (strcmp("W3m0m1", field) == 0)
    field_ind = 75;
    if (strcmp("W3m1m1", field) == 0)
    field_ind = 76;
    if (strcmp("W3m2m1", field) == 0)
    field_ind = 77;
    if (strcmp("W3p3m2", field) == 0)
    field_ind = 78;
    if (strcmp("W3p2m2", field) == 0)
    field_ind = 79;
    if (strcmp("W3p1m2", field) == 0)
    field_ind = 80;
    if (strcmp("W3m0m2", field) == 0)
    field_ind = 81;
    if (strcmp("W3m1m2", field) == 0)
    field_ind = 82;
    if (strcmp("W3m2m2", field) == 0)
    field_ind = 83;
    if (strcmp("W3p3m3", field) == 0)
    field_ind = 84;
    if (strcmp("W3p2m3", field) == 0)
    field_ind = 85;
    if (strcmp("W3p1m3", field) == 0)
    field_ind = 86;
    if (strcmp("W3m0m3", field) == 0)
    field_ind = 87;
    if (strcmp("W3m1m3", field) == 0)
    field_ind = 88;
    if (strcmp("W3m2m3", field) == 0)
    field_ind = 89;
    if (strcmp("W3m3m3", field) == 0)
    field_ind = 90;
    if (strcmp("W3m3m2", field) == 0)
    field_ind = 91;
    if (strcmp("W3m3m1", field) == 0)
    field_ind = 92;
    if (strcmp("W3m3m0", field) == 0)
    field_ind = 93;
    if (strcmp("W3m3p1", field) == 0)
    field_ind = 94;
    if (strcmp("W3m3p2", field) == 0)
    field_ind = 95;
    if (strcmp("W3m3p3", field) == 0)
    field_ind = 96;
    
    // W4
    if (strcmp("W4p2m2", field) == 0 || strcmp("W4p1m2", field) == 0)
    field_ind = 97;
    if (strcmp("W4m0m2", field) == 0 || strcmp("W4m1m2", field) == 0)
    field_ind = 98;
    if (strcmp("W4p2m1", field) == 0 || strcmp("W4p1m1", field) == 0)
    field_ind = 99;
    if (strcmp("W4m0m1", field) == 0 || strcmp("W4m1m1", field) == 0)
    field_ind = 100;
    if (strcmp("W4p2m0", field) == 0 || strcmp("W4p1m0", field) == 0)
    field_ind = 101;
    if (strcmp("W4m0m0", field) == 0 || strcmp("W4m1m0", field) == 0)
    field_ind = 102;
    if (strcmp("W4p1p1", field) == 0 || strcmp("W4m0p1", field) == 0)
    field_ind = 103;
    if (strcmp("W4m1p1", field) == 0 || strcmp("W4m1p2", field) == 0 || strcmp("W4m1p3", field) == 0)
    field_ind = 104;
    if (strcmp("W4m2m0", field) == 0 || strcmp("W4m2p1", field) == 0)
    field_ind = 105;
    if (strcmp("W4m2p2", field) == 0 || strcmp("W4m2p3", field) == 0)
    field_ind = 106;
    if (strcmp("W4m3m0", field) == 0 || strcmp("W4m3p1", field) == 0)
    field_ind = 107;
    if (strcmp("W4m3p2", field) == 0 || strcmp("W4m3p3", field) == 0)
    field_ind = 108;
    
    return field_ind;
}

int jackarea3 (char field[]) {
    // A function to determine the jackknife region for the galaxies.
    int field_ind = 0;
    
    // W1
    if (strcmp("W1p4p3", field) == 0 || strcmp("W1p3p3", field) == 0 || strcmp("W1p2p3", field) == 0 || strcmp("W1p1p3", field) == 0 || strcmp("W1m0p3", field) == 0 || strcmp("W1m1p3", field) == 0 || strcmp("W1m2p3", field) == 0 || strcmp("W1m3p3", field) == 0)
    field_ind = 0;
    if (strcmp("W1p4p2", field) == 0 || strcmp("W1p3p2", field) == 0 || strcmp("W1p2p2", field) == 0 || strcmp("W1p1p2", field) == 0 || strcmp("W1m0p2", field) == 0 || strcmp("W1m1p2", field) == 0 || strcmp("W1m2p2", field) == 0 || strcmp("W1m3p2", field) == 0)
    field_ind = 1;
    if (strcmp("W1p4p1", field) == 0 || strcmp("W1p3p1", field) == 0 || strcmp("W1p2p1", field) == 0 || strcmp("W1p1p1", field) == 0 || strcmp("W1m0p1", field) == 0 || strcmp("W1m1p1", field) == 0 || strcmp("W1m2p1", field) == 0 || strcmp("W1m3p1", field) == 0)
    field_ind = 2;
    if (strcmp("W1p4m0", field) == 0 || strcmp("W1p3m0", field) == 0 || strcmp("W1p2m0", field) == 0 || strcmp("W1p1m0", field) == 0 || strcmp("W1m0m0", field) == 0 || strcmp("W1m1m0", field) == 0 || strcmp("W1m2m0", field) == 0 || strcmp("W1m3m0", field) == 0)
    field_ind = 3;
    if (strcmp("W1p4m1", field) == 0 || strcmp("W1p3m1", field) == 0 || strcmp("W1p2m1", field) == 0 || strcmp("W1p1m1", field) == 0 || strcmp("W1m0m1", field) == 0 || strcmp("W1m1m1", field) == 0 || strcmp("W1m2m1", field) == 0 || strcmp("W1m3m1", field) == 0)
    field_ind = 4;
    if (strcmp("W1p4m2", field) == 0 || strcmp("W1p3m2", field) == 0 || strcmp("W1p2m2", field) == 0 || strcmp("W1p1m2", field) == 0 || strcmp("W1m0m2", field) == 0 || strcmp("W1m1m2", field) == 0 || strcmp("W1m2m2", field) == 0 || strcmp("W1m3m2", field) == 0)
    field_ind = 5;
    if (strcmp("W1p4m3", field) == 0 || strcmp("W1p3m3", field) == 0 || strcmp("W1p2m3", field) == 0 || strcmp("W1p1m3", field) == 0 || strcmp("W1m0m3", field) == 0 || strcmp("W1m1m3", field) == 0 || strcmp("W1m2m3", field) == 0 || strcmp("W1m3m3", field) == 0)
    field_ind = 6;
    if (strcmp("W1p4m4", field) == 0 || strcmp("W1p3m4", field) == 0 || strcmp("W1p2m4", field) == 0 || strcmp("W1p1m4", field) == 0 || strcmp("W1m0m4", field) == 0 || strcmp("W1m1m4", field) == 0 || strcmp("W1m2m4", field) == 0 || strcmp("W1m3m4", field) == 0)
    field_ind = 7;
    if (strcmp("W1m4p3", field) == 0 || strcmp("W1m4p2", field) == 0 || strcmp("W1m4p1", field) == 0 || strcmp("W1m4m0", field) == 0 || strcmp("W1m4m1", field) == 0 || strcmp("W1m4m2", field) == 0 || strcmp("W1m4m3", field) == 0 || strcmp("W1m4m4", field) == 0)
    field_ind = 8;
    
    // W2
    if (strcmp("W2p3p2", field) == 0 || strcmp("W2p3p1", field) == 0 || strcmp("W2p3m0", field) == 0 || strcmp("W2p3m1", field) == 0 || strcmp("W2p2p2", field) == 0 || strcmp("W2p2p1", field) == 0 || strcmp("W2p2m0", field) == 0 || strcmp("W2p2m1", field) == 0)
    field_ind = 9;
    if (strcmp("W2p1p2", field) == 0 || strcmp("W2p1p1", field) == 0 || strcmp("W2p1m0", field) == 0 || strcmp("W2p1m1", field) == 0 || strcmp("W2m0p2", field) == 0 || strcmp("W2m0p1", field) == 0 || strcmp("W2m0m0", field) == 0 || strcmp("W2m0m1", field) == 0)
    field_ind = 10;
    if (strcmp("W2m1p2", field) == 0 || strcmp("W2m1p1", field) == 0 || strcmp("W2m1m0", field) == 0 || strcmp("W2m1m1", field) == 0 || strcmp("W2p3p3", field) == 0 || strcmp("W2p2p3", field) == 0 || strcmp("W2p1p3", field) == 0 || strcmp("W2m0p3", field) == 0 || strcmp("W2m1p3", field) == 0)
    field_ind = 11;
    
    // W3
    if (strcmp("W3p3p3", field) == 0 || strcmp("W3p2p3", field) == 0 || strcmp("W3p1p3", field) == 0 || strcmp("W3m0p3", field) == 0)
    field_ind = 12;
    if (strcmp("W3m1p3", field) == 0 || strcmp("W3m2p3", field) == 0 || strcmp("W3p3p2", field) == 0 || strcmp("W3p2p2", field) == 0)
    field_ind = 13;
    if (strcmp("W3p1p2", field) == 0 || strcmp("W3m0p2", field) == 0 || strcmp("W3m1p2", field) == 0 || strcmp("W3m2p2", field) == 0)
    field_ind = 14;
    if (strcmp("W3p3p1", field) == 0 || strcmp("W3p2p1", field) == 0 || strcmp("W3p1p1", field) == 0 || strcmp("W3m0p1", field) == 0)
    field_ind = 15;
    if (strcmp("W3m1p1", field) == 0 || strcmp("W3m2p1", field) == 0 || strcmp("W3p3m0", field) == 0 || strcmp("W3p2m0", field) == 0)
    field_ind = 16;
    if (strcmp("W3p1m0", field) == 0 || strcmp("W3m0m0", field) == 0 || strcmp("W3m1m0", field) == 0 || strcmp("W3m2m0", field) == 0)
    field_ind = 17;
    if (strcmp("W3p3m1", field) == 0 || strcmp("W3p2m1", field) == 0 || strcmp("W3p1m1", field) == 0 || strcmp("W3m0m1", field) == 0)
    field_ind = 18;
    if (strcmp("W3m1m1", field) == 0 || strcmp("W3m2m1", field) == 0 || strcmp("W3p3m2", field) == 0 || strcmp("W3p2m2", field) == 0)
    field_ind = 19;
    if (strcmp("W3p1m2", field) == 0 || strcmp("W3m0m2", field) == 0 || strcmp("W3m1m2", field) == 0 || strcmp("W3m2m2", field) == 0)
    field_ind = 20;
    if (strcmp("W3p3m3", field) == 0 || strcmp("W3p2m3", field) == 0 || strcmp("W3p1m3", field) == 0 || strcmp("W3m0m3", field) == 0)
    field_ind = 21;
    if (strcmp("W3m1m3", field) == 0 || strcmp("W3m2m3", field) == 0 || strcmp("W3m3m3", field) == 0 || strcmp("W3m3m2", field) == 0)
    field_ind = 22;
    if (strcmp("W3m3m1", field) == 0 || strcmp("W3m3m0", field) == 0 || strcmp("W3m3p1", field) == 0 || strcmp("W3m3p2", field) == 0 || strcmp("W3m3p3", field) == 0)
    field_ind = 23;
    
    // W4
    if (strcmp("W4p2m2", field) == 0 || strcmp("W4p1m2", field) == 0 || strcmp("W4m0m2", field) == 0 || strcmp("W4m1m2", field) == 0 || strcmp("W4p2m1", field) == 0 || strcmp("W4p1m1", field) == 0 || strcmp("W4m0m1", field) == 0 || strcmp("W4m1m1", field) == 0)
    field_ind = 24;
    if (strcmp("W4p2m0", field) == 0 || strcmp("W4p1m0", field) == 0 || strcmp("W4m0m0", field) == 0 || strcmp("W4m1m0", field) == 0 || strcmp("W4p1p1", field) == 0 || strcmp("W4m0p1", field) == 0 || strcmp("W4m1p1", field) == 0 || strcmp("W4m1p2", field) == 0 || strcmp("W4m1p3", field) == 0)
    field_ind = 25;
    if (strcmp("W4m2m0", field) == 0 || strcmp("W4m2p1", field) == 0 || strcmp("W4m2p2", field) == 0 || strcmp("W4m2p3", field) == 0 || strcmp("W4m3m0", field) == 0 || strcmp("W4m3p1", field) == 0 || strcmp("W4m3p2", field) == 0 || strcmp("W4m3p3", field) == 0)
    field_ind = 26;
    
    return field_ind;
}


void skip_line (int spot, FILE *readfile, int line, int skipped_line)
{
    // A function to skip certain lines for a file.
    do {
        spot = fgetc (readfile);
        if (spot == '\n')
            line++;
        if (line == skipped_line)
            break;
    } while (spot != EOF);
    
    return;
}

int min(double arr[], int size, double subtract)
{
    // A function to find the minimum value of the difference of the array elements to the subtraction.
    double minimum = fabs(arr[0] - subtract);
    int j = 0, count = 0;
    for (j = 0; j < size; j++) {
        if (minimum > fabs(arr[j] - subtract)) {
            minimum = fabs(arr[j] - subtract);
            count = j;
        }
    }
    return count;
}

int max(double arr[], int size, double subtract)
{
    // A function to find the maximum value of the difference of the array elements to the subtraction.
    double maximum = fabs(arr[0] - subtract);
    int j = 0, count = 0;
    for (j = 0; j < size; j++) {
        if (maximum < fabs(arr[j] - subtract)) {
            maximum = fabs(arr[j] - subtract);
            count = j;
        }
    }
    return count;
}

void min_matrix(int row, int column, double **arr, int *row_ind, int *col_ind, double subtract)
{
    // A function to find the minimum value of the difference of the array elements to the subtraction.
    double minimum = fabs(arr[0][0] - subtract);
    int j = 0, k = 0;
    for (j = 0; j < row; j++) {
        for (k = 0; k < column; k++) {
            if (minimum > fabs(arr[j][k] - subtract)) {
                minimum = fabs(arr[j][k] - subtract);
                *row_ind = j;
                *col_ind = k;
            }
        }
    }
}

int near_dec(int column, double *dec_temp, double gal_dec) {
    // A function to find the nearest dec value where the random locates.
    double minimum = fabs(dec_temp[0] - gal_dec);
    int dec_ind = 0;
    int i = 0;
    for (i = 0; i < column; i++) {
        if (minimum > fabs(dec_temp[i] - gal_dec)) {
            minimum = fabs(dec_temp[i] - gal_dec);
            dec_ind = i;
        }
    }
    return dec_ind;
}

int near_ra(int column, double *ra_temp, double gal_ra) {
    // A function to find the nearest ra value where the random locates.
    double minimum = fabs(ra_temp[0] - gal_ra);
    int ra_ind = 0;
    int j = 0;
    for (j = 0; j < column; j++) {
        if (minimum > fabs(ra_temp[j] - gal_ra)) {
            minimum = fabs(ra_temp[j] - gal_ra);
            ra_ind = j;
        }
    }
    return ra_ind;
}

double min_array(double arr[], int size) {
    // A function to find the smallest array value.
    double min_d = arr[0];
    int i = 0;
    for (i = 1; i < size; i++) {
        if (arr[i] < min_d)
            min_d = arr[i];
    }
    return min_d;
}

double max_array(double arr[], int size) {
    // A function to find the largest array value.
    double max_d = arr[0];
    int i = 0;
    for (i = 1; i < size; i++) {
        if (arr[i] > max_d)
            max_d = arr[i];
    }
    return max_d;
}

double min_value(double arr[], int size)
{
    // A function to find the minimum value of the array
    double minimum = arr[0];
    int j = 0, count = 0;
    for (j = 0; j < size; j++) {
        if (minimum > arr[j]) {
            minimum = arr[j];
        }
    }
    return minimum;
}

double max_value(double arr[], int size)
{
    // A function to find the maximum value of the array
    double maximum = arr[0];
    int j = 0, count = 0;
    for (j = 0; j < size; j++) {
        if (maximum < arr[j]) {
            maximum = arr[j];
        }
    }
    return maximum;
}

int max_int(int arr[], int size)
{
    // A function to find the maximum integer of the array
    int maximum = arr[0];
    int j = 0, count = 0;
    for (j = 0; j < size; j++) {
        if (maximum < arr[j]) {
            maximum = arr[j];
        }
    }
    return maximum;
}

double arr_mean(double arr[], int size)
{
    // A function to compute the mean of the array elements.
    double avg = 0, total = 0;
    int i = 0;
    for (i = 0; i < size; i++) {
        total += arr[i];
    }
    avg = total / size;
    return avg;
}

double arr_var(double arr[], double avg, int size)
{
    // A function to compute the sample variance of the array elements.
    double var = 0, total = 0;
    int i = 0;
    for (i = 0; i < size; i++) {
        total += pow((arr[i] - avg), 2);
    }
    var = total / (size - 1);
    return var;
}

double max_compare(double a, double b) {
    // A function to find the MAX{a, b}.
    if (a >= b) return a;
    else    return b;
}

double dist_periodic(double x1, double y1, double x2, double y2, double L) {
    // Calculate the distance between (x1, y1) and (x2, y2) with periodic boundary condition.
    return hypot(remainder(x1-x2,L),remainder(y1-y2,L));
}

int ragrid_number(double r_search, double delta, double rasize) {
    // A function to find the number of grids along RA, with celestial geometry accounted for. Input of r_search should be in arcsecs, delta should be in deg, rasize should be in arcsecs.
    double delta_radian = delta * M_PI / 180.;
    double r_radian = r_search * M_PI / (180 * 3600.);
    double size = 0;
    int ra_num = 0;
    
    size = acos((cos(r_radian) - pow(sin(delta_radian), 2)) / pow(cos(delta_radian), 2)) * (3600 * 180) / M_PI;
    ra_num = floor(rasize / size);
    return ra_num;
}

int compare ( const void *pa, const void *pb ) {
    // Sorting of the n * 3 array. 1st col in ascending order, if 2nd col (DEC) is degenerate, then sort 1st (RA), if still degenerate, then 3rd.
    const int *a = *(const int **)pa;
    const int *b = *(const int **)pb;
    if (a[1] < b[1]) return -1;
    if (a[1] > b[1]) return 1;
    if (a[1] == b[1]) {
        if (a[0] < b[0]) return -1;
        if (a[0] > b[0]) return 1;
        if (a[0] == b[0]) {
            if (a[2] < b[2]) return -1;
            if (a[2] > b[2]) return 1;
            else return 0;
        }
    }
}

int compare2 (const void * a, const void * b) {
    // Sorting an array.
    return ( *(int*)a - *(int*)b );
}

int cmp(const void *x, const void *y) {
    // Sorting a double array.
    double xx = *(double*)x, yy = *(double*)y;
    if (xx < yy) return -1;
    if (xx > yy) return  1;
    return 0;
}

int grid2consider (int **gridscon, int ra_grid, int ra_mesh, int dec_grid, int dec_mesh) {
    // A function to store in arrays which pixels to consider.
    if (ra_grid == 0 && dec_grid == 0) {
        gridscon[0][0] = 0;
        gridscon[0][1] = 0;
        gridscon[1][0] = 1;
        gridscon[1][1] = 0;
        gridscon[2][0] = 0;
        gridscon[2][1] = 1;
        gridscon[3][0] = 1;
        gridscon[3][1] = 1;
    }
    else if (ra_grid == ra_mesh - 1 && dec_grid == 0) {
        gridscon[0][0] = ra_mesh - 2;
        gridscon[0][1] = 0;
        gridscon[1][0] = ra_mesh - 1;
        gridscon[1][1] = 0;
        gridscon[2][0] = ra_mesh - 2;
        gridscon[2][1] = 1;
        gridscon[3][0] = ra_mesh - 1;
        gridscon[3][1] = 1;
    }
    else if (ra_grid == 0 && dec_grid == dec_mesh - 1) {
        gridscon[0][0] = 0;
        gridscon[0][1] = dec_mesh - 2;
        gridscon[1][0] = 1;
        gridscon[1][1] = dec_mesh - 2;
        gridscon[2][0] = 0;
        gridscon[2][1] = dec_mesh - 1;
        gridscon[3][0] = 1;
        gridscon[3][1] = dec_mesh - 1;
    }
    else if (ra_grid == ra_mesh - 1 && dec_grid == dec_mesh - 1) {
        gridscon[0][0] = ra_mesh - 2;
        gridscon[0][1] = dec_mesh - 2;
        gridscon[1][0] = ra_mesh - 1;
        gridscon[1][1] = dec_mesh - 2;
        gridscon[2][0] = ra_mesh - 2;
        gridscon[2][1] = dec_mesh -1;
        gridscon[3][0] = ra_mesh - 1;
        gridscon[3][1] = dec_mesh - 1;
    }
    // Margins, not the corners.
    else if (ra_grid == 0 && dec_grid != 0 && dec_grid != dec_mesh - 1) {
        gridscon[0][0] = 0;
        gridscon[0][1] = dec_grid - 1;
        gridscon[1][0] = 1;
        gridscon[1][1] = dec_grid - 1;
        gridscon[2][0] = 0;
        gridscon[2][1] = dec_grid;
        gridscon[3][0] = 1;
        gridscon[3][1] = dec_grid;
        gridscon[4][0] = 0;
        gridscon[4][1] = dec_grid + 1;
        gridscon[5][0] = 1;
        gridscon[5][1] = dec_grid + 1;
    }
    else if (ra_grid == ra_mesh - 1 && dec_grid != 0 && dec_grid != dec_mesh - 1) {
        gridscon[0][0] = ra_mesh - 2;
        gridscon[0][1] = dec_grid - 1;
        gridscon[1][0] = ra_mesh - 1;
        gridscon[1][1] = dec_grid - 1;
        gridscon[2][0] = ra_mesh - 2;
        gridscon[2][1] = dec_grid;
        gridscon[3][0] = ra_mesh - 1;
        gridscon[3][1] = dec_grid;
        gridscon[4][0] = ra_mesh - 2;
        gridscon[4][1] = dec_grid + 1;
        gridscon[5][0] = ra_mesh - 1;
        gridscon[5][1] = dec_grid + 1;
    }
    else if (ra_grid != 0 && ra_grid != ra_mesh - 1 && dec_grid == 0) {
        gridscon[0][0] = ra_grid - 1;
        gridscon[0][1] = 0;
        gridscon[1][0] = ra_grid;
        gridscon[1][1] = 0;
        gridscon[2][0] = ra_grid + 1;
        gridscon[2][1] = 0;
        gridscon[3][0] = ra_grid - 1;
        gridscon[3][1] = 1;
        gridscon[4][0] = ra_grid;
        gridscon[4][1] = 1;
        gridscon[5][0] = ra_grid + 1;
        gridscon[5][1] = 1;
    }
    else if (ra_grid != 0 && ra_grid != ra_mesh - 1 && dec_grid == dec_mesh - 1) {
        gridscon[0][0] = ra_grid - 1;
        gridscon[0][1] = dec_mesh - 2;
        gridscon[1][0] = ra_grid;
        gridscon[1][1] = dec_mesh - 2;
        gridscon[2][0] = ra_grid + 1;
        gridscon[2][1] = dec_mesh - 2;
        gridscon[3][0] = ra_grid - 1;
        gridscon[3][1] = dec_mesh - 1;
        gridscon[4][0] = ra_grid;
        gridscon[4][1] = dec_mesh - 1;
        gridscon[5][0] = ra_grid + 1;
        gridscon[5][1] = dec_mesh - 1;
    }
    // General case.
    else if (ra_grid != 0 && ra_grid != ra_mesh - 1 && dec_grid != 0 && dec_grid != dec_mesh - 1) {
        gridscon[0][0] = ra_grid - 1;
        gridscon[0][1] = dec_grid - 1;
        gridscon[1][0] = ra_grid;
        gridscon[1][1] = dec_grid - 1;
        gridscon[2][0] = ra_grid + 1;
        gridscon[2][1] = dec_grid - 1;
        gridscon[3][0] = ra_grid - 1;
        gridscon[3][1] = dec_grid;
        gridscon[4][0] = ra_grid;
        gridscon[4][1] = dec_grid;
        gridscon[5][0] = ra_grid + 1;
        gridscon[5][1] = dec_grid;
        gridscon[6][0] = ra_grid - 1;
        gridscon[6][1] = dec_grid + 1;
        gridscon[7][0] = ra_grid;
        gridscon[7][1] = dec_grid + 1;
        gridscon[8][0] = ra_grid + 1;
        gridscon[8][1] = dec_grid + 1;
    }

    return 0;
}

double ang_sep (double ra_cen, double dec_cen, double ra, double dec) {
    double theta = 0;
    double ra_diff = (ra - ra_cen) * M_PI / 180.;
    double dec1 = dec_cen * M_PI / 180.;
    double dec2 = dec * M_PI / 180.;
    theta = acos(cos(ra_diff) * cos(dec1) * cos(dec2) + sin(dec1) * sin(dec2)) * 180. * 3600 / M_PI;
    // Theta in arcsecs.
    return theta;
}

double fRand(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

double ang_sep_rad (double ra_cen, double dec_cen, double ra, double dec) {
    double theta = 0;
    double ra_diff = (ra - ra_cen) * M_PI / 180.;
    double dec1 = dec_cen * M_PI / 180.;
    double dec2 = dec * M_PI / 180.;
    theta = acos(cos(ra_diff) * cos(dec1) * cos(dec2) + sin(dec1) * sin(dec2));
    // Theta in radian.
    return theta;
}

double ori_ang (double a1, double d1, double a2, double d2, double ang_sep) {
    // The orientation of the source galaxy w.r.t. the positive declination.
    // Input a1, a2, d1, d2 in deg. ang_sep in arcsecs.
    // 1 and 2 are lens and source galaxies respectively.
    double d2r = M_PI / 180.;
    double r2d = 180. / M_PI;
    double ra1 = a1 * d2r;
    double ra2 = a2 * d2r;
    double dec1 = d1 * d2r;
    double dec2 = d2 * d2r;
    double a_sep = ang_sep * d2r / 3600;
    double o_ang = 0;
    double mod = fabs(cos(dec1) * sin(ra2 - ra1) / sin(a_sep));
    
    if (mod > 1) mod = 1;
    
    o_ang = acos(mod);
    
    // Angle between the great circle and the declination of the source. [0, pi] for the spin-2 field.
    if (a1 <= a2 && d1 < d2)    return o_ang;   // Source in 1st quadrant.
    if (a1 > a2 && d1 <= d2)    return M_PI - o_ang;    // 2nd quadrant.
    if (a1 >= a2 && d1 > d2)    return o_ang;   // 3rd quadrant.
    if (a1 < a2 && d1 >= d2)    return M_PI - o_ang;    // 4th quadrant.
}

int lenslumclass (char* lenscolor, char* lum_class) {
    // Determine the color and luminosity class of the lens galaxy.
    if (strcmp(lenscolor, "blu") == 0) {
        if (strcmp(lum_class, "L1") == 0) return 0;
        if (strcmp(lum_class, "L2") == 0) return 1;
        if (strcmp(lum_class, "L3") == 0) return 2;
        if (strcmp(lum_class, "L4") == 0) return 3;
        if (strcmp(lum_class, "L5") == 0) return 4;
        if (strcmp(lum_class, "L6") == 0) return 5;
        if (strcmp(lum_class, "L7") == 0) return 6;
        if (strcmp(lum_class, "L8") == 0) return 7;
    }
    if (strcmp(lenscolor, "red") == 0) {
        if (strcmp(lum_class, "L1") == 0) return 8;
        if (strcmp(lum_class, "L2") == 0) return 9;
        if (strcmp(lum_class, "L3") == 0) return 10;
        if (strcmp(lum_class, "L4") == 0) return 11;
        if (strcmp(lum_class, "L5") == 0) return 12;
        if (strcmp(lum_class, "L6") == 0) return 13;
        if (strcmp(lum_class, "L7") == 0) return 14;
        if (strcmp(lum_class, "L8") == 0) return 15;
    }
}

int lenslumclass_wider (char* lenscolor, char* lum_class) {
    // Determine the color and luminosity class of the lens galaxy.
    if (strcmp(lenscolor, "blu") == 0) {
        if (strcmp(lum_class, "L1") == 0) return 0;
        if (strcmp(lum_class, "L2") == 0) return 1;
        if (strcmp(lum_class, "L3") == 0) return 2;
        if (strcmp(lum_class, "L4") == 0) return 3;
    }
    if (strcmp(lenscolor, "red") == 0) {
        if (strcmp(lum_class, "L1") == 0) return 4;
        if (strcmp(lum_class, "L2") == 0) return 5;
        if (strcmp(lum_class, "L3") == 0) return 6;
        if (strcmp(lum_class, "L4") == 0) return 7;
    }
}


int lensmulticlass (char* lenscolor, int color_bin, char* lum_class, int lum_bin, char* size_class, int size_bin, char* mass_class, int mass_bin) {
    
    int total_bin = color_bin * lum_bin * size_bin * mass_bin;
    int color_index = 0, lum_index = 0, size_index = 0, mass_index = 0, total_index = 0, i = 0, j = 0, k = 0, l = 0;
    
    // Determine the color, luminosity, size and mass classes of the lens galaxy.
    if (strcmp(lenscolor, "blu") == 0)  color_index = 0;
    if (strcmp(lenscolor, "red") == 0)  color_index = 1;
    
    if (strcmp(lum_class, "L1") == 0)  lum_index = 0;
    if (strcmp(lum_class, "L2") == 0)  lum_index = 1;
    if (strcmp(lum_class, "L3") == 0)  lum_index = 2;
    if (strcmp(lum_class, "L4") == 0)  lum_index = 3;
    if (strcmp(lum_class, "L5") == 0)  lum_index = 4;
    if (strcmp(lum_class, "L6") == 0)  lum_index = 5;
    if (strcmp(lum_class, "L7") == 0)  lum_index = 6;
    
    if (strcmp(size_class, "S1") == 0)  size_index = 0;
    if (strcmp(size_class, "S2") == 0)  size_index = 1;
    if (strcmp(size_class, "S3") == 0)  size_index = 2;
    if (strcmp(size_class, "S4") == 0)  size_index = 3;
    if (strcmp(size_class, "S5") == 0)  size_index = 4;
    if (strcmp(size_class, "S6") == 0)  size_index = 5;
    
    if (strcmp(mass_class, "M1") == 0)  mass_index = 0;
    if (strcmp(mass_class, "M2") == 0)  mass_index = 1;
    if (strcmp(mass_class, "M3") == 0)  mass_index = 2;
    if (strcmp(mass_class, "M4") == 0)  mass_index = 3;
    if (strcmp(mass_class, "M5") == 0)  mass_index = 4;
    if (strcmp(mass_class, "M6") == 0)  mass_index = 5;
    
    total_index = color_index * (lum_bin * size_bin * mass_bin) + lum_index * (size_bin * mass_bin) + size_index * mass_bin + mass_index;
    
    if (total_index >= total_bin) {
        printf("Class index larger than or equal to total number of bins!\n");
        exit(1);
    }
    
    return total_index;
}


void lens_filename (int lenslumclass, char name[]) {
    // Given lenslumclass, return the name of the file. It corresponds to the function lenslumclass.
    if (lenslumclass == 0) strcpy (name, "blu_L1.dat");
    if (lenslumclass == 1) strcpy (name, "blu_L2.dat");
    if (lenslumclass == 2) strcpy (name, "blu_L3.dat");
    if (lenslumclass == 3) strcpy (name, "blu_L4.dat");
    if (lenslumclass == 4) strcpy (name, "blu_L5.dat");
    if (lenslumclass == 5) strcpy (name, "blu_L6.dat");
    if (lenslumclass == 6) strcpy (name, "blu_L7.dat");
    if (lenslumclass == 7) strcpy (name, "blu_L8.dat");
    if (lenslumclass == 8) strcpy (name, "red_L1.dat");
    if (lenslumclass == 9) strcpy (name, "red_L2.dat");
    if (lenslumclass == 10) strcpy (name, "red_L3.dat");
    if (lenslumclass == 11) strcpy (name, "red_L4.dat");
    if (lenslumclass == 12) strcpy (name, "red_L5.dat");
    if (lenslumclass == 13) strcpy (name, "red_L6.dat");
    if (lenslumclass == 14) strcpy (name, "red_L7.dat");
    if (lenslumclass == 15) strcpy (name, "red_L8.dat");
}

void lens_filename2 (int lenslumclass, char name[]) {
    // Given lenslumclass, return the name of the file. It corresponds to the function lenslumclass.
    if (lenslumclass == 0) strcpy (name, "blu_L1_");
    if (lenslumclass == 1) strcpy (name, "blu_L2_");
    if (lenslumclass == 2) strcpy (name, "blu_L3_");
    if (lenslumclass == 3) strcpy (name, "blu_L4_");
    if (lenslumclass == 4) strcpy (name, "blu_L5_");
    if (lenslumclass == 5) strcpy (name, "blu_L6_");
    if (lenslumclass == 6) strcpy (name, "blu_L7_");
    if (lenslumclass == 7) strcpy (name, "blu_L8_");
    if (lenslumclass == 8) strcpy (name, "red_L1_");
    if (lenslumclass == 9) strcpy (name, "red_L2_");
    if (lenslumclass == 10) strcpy (name, "red_L3_");
    if (lenslumclass == 11) strcpy (name, "red_L4_");
    if (lenslumclass == 12) strcpy (name, "red_L5_");
    if (lenslumclass == 13) strcpy (name, "red_L6_");
    if (lenslumclass == 14) strcpy (name, "red_L7_");
    if (lenslumclass == 15) strcpy (name, "red_L8_");
}

void lens_filename_assembias (int lenslumclass, char name[]) {
    // Given lenslumclass, return the name of the file. It corresponds to the function lenslumclass.
    if (lenslumclass == 0) strcpy (name, "blu_L1_");
    if (lenslumclass == 1) strcpy (name, "blu_L2_");
    if (lenslumclass == 2) strcpy (name, "red_L1_");
    if (lenslumclass == 3) strcpy (name, "red_L2_");
}

void lens_filenamemultifx (int lenslumclass, char name[], int color_bin, int lum_bin, int size_bin, int mass_bin) {
    
    // Given lenslumclass, return the name of the file. It corresponds to the function lenslumclass.
    
    int color_rem = 0, lum_rem = 0, size_rem = 0, sub1 = 0, sub2 = 0, sub3 = 0;
    
    color_rem = lenslumclass / (lum_bin * size_bin * mass_bin);
    
    if (color_rem == 0) strcpy (name, "blu_");
    if (color_rem == 1) strcpy (name, "red_");
    
    sub1 = lenslumclass % (lum_bin * size_bin * mass_bin);
    lum_rem = sub1 / (size_bin * mass_bin);
    
    if (lum_rem == 0) strcat (name, "L1_");
    if (lum_rem == 1) strcat (name, "L2_");
    if (lum_rem == 2) strcat (name, "L3_");
    if (lum_rem == 3) strcat (name, "L4_");
    if (lum_rem == 4) strcat (name, "L5_");
    if (lum_rem == 5) strcat (name, "L6_");
    if (lum_rem == 6) strcat (name, "L7_");
    
    sub2 = sub1 % (size_bin * mass_bin);
    size_rem = sub2 / mass_bin;
    
    if (size_rem == 0) strcat (name, "S1_");
    if (size_rem == 1) strcat (name, "S2_");
    if (size_rem == 2) strcat (name, "S3_");
    if (size_rem == 3) strcat (name, "S4_");
    if (size_rem == 4) strcat (name, "S5_");
    if (size_rem == 5) strcat (name, "S6_");
    
    sub3 = sub2 % mass_bin;
    
    if (sub3 == 0) strcat (name, "M1_");
    if (sub3 == 1) strcat (name, "M2_");
    if (sub3 == 2) strcat (name, "M3_");
    if (sub3 == 3) strcat (name, "M4_");
    if (sub3 == 4) strcat (name, "M5_");
    if (sub3 == 5) strcat (name, "M6_");
}

void lens_filename2_wider (int lenslumclass_wider, char name[]) {
    // Given lenslumclass, return the name of the file. It corresponds to the function lenslumclass.
    if (lenslumclass_wider == 0) strcpy (name, "blu_L1_");
    if (lenslumclass_wider == 1) strcpy (name, "blu_L2_");
    if (lenslumclass_wider == 2) strcpy (name, "blu_L3_");
    if (lenslumclass_wider == 3) strcpy (name, "blu_L4_");
    if (lenslumclass_wider == 4) strcpy (name, "red_L1_");
    if (lenslumclass_wider == 5) strcpy (name, "red_L2_");
    if (lenslumclass_wider == 6) strcpy (name, "red_L3_");
    if (lenslumclass_wider == 7) strcpy (name, "red_L4_");
}

