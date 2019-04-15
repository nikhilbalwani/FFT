#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.1415927
#endif

typedef struct complex_num_struct {
    double real;
    double imag;
} complex_num;

complex_num add(complex_num num1, complex_num num2);
complex_num subtract(complex_num num1, complex_num num2);
complex_num multiply(complex_num num1, complex_num num2);
complex_num divide(complex_num num1, complex_num num2);
complex_num power(complex_num num, double n);
complex_num complex_exp(double theta);
double magnitude(complex_num num);
double phase(complex_num num);
void display_complex_matrix(complex_num** matrix, int height, int width);
void fft2(complex_num** input, complex_num** output, int height, int width);
void fft_driver(complex_num input[], complex_num output[], int n, int step);
void fft(complex_num input[], complex_num * output[], int n);
void display_complex_vec(complex_num * vec, int n);
complex_num conjugate(complex_num num);
complex_num ** matrix;
complex_num ** temp;

int main() {

    complex_num one;
    one.real = 1;
    one.imag = 0;

    complex_num height, width;

    height.real = 4;
    width.real = 4;
    height.imag = 0;
    width.imag = 0;

    complex_num zero;
    zero.real = 0;
    zero.imag = 0;

    complex_num input1[] = {one, one, one, one};
    complex_num input2[] = {one, one, one, zero};
    complex_num input3[] = {one, one, zero, zero};
    complex_num input4[] = {one, zero, zero, zero};

    complex_num* input[] = {input1, input2, input3, input4};
    complex_num* output[] = {input1, input2, input3, input4};
    complex_num* conj[] = {input1, input2, input3, input4};

    fft2(input, output, 4, 4);
    display_complex_matrix(output, 4, 4);

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            conj[i][j] = conjugate(output[i][j]);
            output[i][j] = conj[i][j];
        }
    }

    fft2(conj, output, 4, 4);

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            conj[i][j] = conjugate(divide(output[i][j], (multiply(height, width) ) ) );
        }
    }

    display_complex_matrix(conj, 4, 4);
}

void fft2(complex_num** input, complex_num** output, int height, int width) {

    matrix = (complex_num **) malloc(height * sizeof(complex_num*));
    temp = (complex_num **) malloc(width * sizeof(complex_num*));

    for (int i = 0; i < height; ++i) {
        fft(input[i], &temp[i], width);
    }

    for (int i = 0; i < height; ++i) {

        matrix[i] = (complex_num*) malloc(width * sizeof(complex_num));

        for (int j = 0; j < width; ++j) {
            matrix[i][j] = temp[j][i];
            printf("%lf", temp[j][i].real);

        }
    }

    for (int i = 0; i < height; ++i) {
        free(temp[i]);
        fft(matrix[i], &temp[i], width);
        free(matrix[i]);
    }

    free(matrix);

    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            output[i][j] = temp[j][i];
        }
    }
    for (int i = 0; i < height; ++i)
      free(temp[i]);
    free(temp);
}

void display_complex_vec(complex_num vec[], int n) {

    printf("[");

    for (int i = 0; i < n; ++i) {
        printf("(%g, %g), ", vec[i].real, vec[i].imag);
    }

    printf("]\n");
}

void display_complex_matrix(complex_num** matrix, int height, int width) {


    for (int i = 0; i < height; ++i) {
        printf("[");

        for (int j = 0; j < width; ++j) {
            printf("(%g, %g), ", matrix[i][j].real, matrix[i][j].imag);
        }

        printf("]\n");
    }
}

void fft_driver(complex_num input[], complex_num output[], int n, int step) {
    complex_num diff, cexp;

    if (step < n) {

        fft_driver(output, input, n, step * 2);
        fft_driver(output + step, input + step, n, step * 2);

        for (int i = 0; i < n; i += 2 * step) {
            cexp = complex_exp(-M_PI * (double) i / (double) n);

            diff = multiply(complex_exp(-M_PI * (double) i / (double) n), output[i + step]);

            input[i / 2] = add(output[i], diff);
            input[(i + n) / 2] = subtract(output[i], diff);
            //printf("(%g, %g) * (%g, %g)\n", input[i / 2].real, input[i / 2].imag, input[(i + n) / 2].real, input[(i + n) / 2].imag);
        }
    }
}

void fft(complex_num* input, complex_num** output, int n) {

  complex_num* dummy = (complex_num *)malloc(n * sizeof(complex_num));
  complex_num* inp = (complex_num *)malloc(n * sizeof(complex_num));

  for (int i = 0; i < n; ++i) {
    dummy[i].real = 0;
    dummy[i].imag = 0;
    inp[i].real = 0;
    inp[i].imag = 0;
  }

  printf("%lf", inp[1].real);
  for (int i = 0; i < n; ++i) {
    dummy[i] = input[i];
    inp[i] = input[i];
  }

  fft_driver(inp, dummy, n, 1);

  free(dummy);
  printf("Hello\n");

  *output = inp;
}

complex_num complex_exp(double theta) {
    complex_num result;
    result.real = cos(theta);
    result.imag = sin(theta);

    return result;
}

double magnitude(complex_num num) {
    return sqrt(pow(num.real, 2) + pow(num.imag, 2));
}

double phase(complex_num num) {
    return atan(num.imag / num.real);
}

complex_num add(complex_num num1, complex_num num2) {
    complex_num result;
    result.real = num1.real + num2.real;
    result.imag = num1.imag + num2.imag;

    return result;
}

complex_num conjugate(complex_num num) {
    complex_num result;
    result.real = num.real;
    result.imag = -num.imag;

    return result;
}

complex_num subtract(complex_num num1, complex_num num2) {
    complex_num result;
    result.real = num1.real - num2.real;
    result.imag = num1.imag - num2.imag;

    return result;
}

complex_num multiply(complex_num num1, complex_num num2) {
    complex_num result;
    result.real = num1.real * num2.real - num1.imag * num2.imag;
    result.imag = num1.real * num2.imag + num1.imag * num2.real;

    return result;
}

complex_num divide(complex_num num1, complex_num num2) {
    complex_num result;

    double magnitude_square = pow(num2.real, 2) + pow(num2.imag, 2);
    num2.imag = -num2.imag;

    result = multiply(num1, num2);
    result.real = result.real / magnitude_square;
    result.imag = result.imag / magnitude_square;

    return result;
}

complex_num power(complex_num num, double n) {
    complex_num result;

    double magnitude = sqrt(pow(num.real, 2) + pow(num.imag, 2));

    result.real = pow(magnitude, n) * cos(n * atan(num.imag / num.real));
    result.imag = pow(magnitude, n) * sin(n * atan(num.imag / num.real));

    return result;
}
