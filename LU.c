#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAXMAT 100

typedef float matrice[MAXMAT][MAXMAT];

// given structure for the square matrix of size nb
struct matrfloat {
  int n;
  matrice m;
};
typedef struct matrfloat mat_t;

// function to fill a square matrix of size nb
mat_t create_A(int nb) {
  mat_t m;
  int i, j, k, l, z, zi;
  m.n = nb;
  for (i = 0, k = nb; i < nb; i++, k--) {
    l = k;
    m.m[i][i] = (float)l;
    l >>= 1;
    for (j = i + 1, z = 1, zi = 0; j < nb; j++) {
      if (z == zi) {
        m.m[i][j] = (float)l;
        m.m[j][i] = (float)l;
        l >>= 1;
        z <<= 1;
        zi = 0;
      } else {
        m.m[i][j] = 0.0;
        m.m[j][i] = 0.0;
        zi++;
      }
    }
  }
  return m;
}

// function to fill a vector of size nb
void create_B(float *B, int nb) {

  int i;
  for (i = 0; i < nb; i++) {
    B[i] = rand() % 10;
  }
}

// function to print a matrix
// takes pointer to mat_m as argument, as well as an int to differentiate between printing the initial matrix A and the combined matrix LU
void print_mat(mat_t *m, int lu) {

  int i, j;

  if (lu == 1) {
    printf("\nLU :\n");
  } else {
    printf("A:\n");
  }

  for (i = 0; i < m->n; i++) {
    for (j = 0; j < m->n; j++) {
      if (m->m[i][j] >= 0) {
        printf(" ");
      }
      if (abs(m->m[i][j]) < 10) {
        printf(" ");
      }
      printf(" %.3f", m->m[i][j]);
    }
    printf("\n");
  }
}

// function to switch two rows in matrix A as well as the corresponding values of vector B in the equation AX=B
// takes pointers to structure mat_t (containing matrix A) and vector B as arguments, as well as int values corresponding to the rows of the matrix to be switched
void switch_rows(mat_t *m, float *B, int row1, int row2) {

  int i;
  float temp;
  for (i = 0; i < m->n; i++) {
    temp = m->m[row1][i];
    m->m[row1][i] = m->m[row2][i];
    m->m[row2][i] = temp;
  }

  temp = B[row1];
  B[row1] = B[row2];
  B[row2] = temp;
}

// function to find the greatest value in a column below a given pivot
// takes pointer to structure mat_t (containing matrix A) and an int representing the row at which and below which we check the column for the greatest absolute value
// returns int representing row at which highest absolute value is found
int find_greatest_in_column(mat_t *m, int row) {

  int i, index, high_value;

  high_value = 0;
  index = row;

  for (i = row; i < m->n; i++) {

    if (abs(m->m[i][row]) > high_value) {
      high_value = m->m[i][row];
      index = i;
    }
  }
  return index;
}

// function to print a vector
// takes pointer to vector B, length of B and int to differentiate between print formatting as arguments
void print_vec(float *B, int len, int x) {
  int i;
  if (x == 1) {
    printf("X:\n");
  } else {
    printf("B:\n");
  }
  for (i = 0; i < len; i++) {
    if (B[i] >= 0) {
      printf(" ");
    }
    if (abs(B[i] < 10)) {
      printf(" ");
    }
    printf(" %.3f\n", B[i]);
  }
}

// function that performs a gaussian pivot
// takes pointer to structure mat_t (containing matrix A) and dimensions of pivot as arguments
void gaussian_pivot(mat_t *m, int sub_matrix_dim) {

  int i, j;
  float pivot, quotient;
  int pivot_coord = m->n - sub_matrix_dim;

  pivot = m->m[pivot_coord][pivot_coord];

  // filling in U matrix
  for (i = pivot_coord; i < m->n - 1; i++) {
    if (pivot == 0) {
      printf("\nAttempting division by zero.\n\nI'm sorry Dave, I can't do that.\n\n");
      exit(0);
    }
    quotient = m->m[i + 1][pivot_coord] / pivot;
    for (j = pivot_coord; j < m->n; j++) {
      m->m[i + 1][j] -= m->m[pivot_coord][j] * quotient;
    }
    // filling in L matrix
    m->m[i + 1][pivot_coord] = quotient;
  }
}

// function to perform gaussian elimination
// takes pointers to structure mat_t (containing matrix A) and vector B as arguments
void gaussian_elimination(mat_t *m, float *B) {

  int i;
  // this int will contain the row to be switched with the pivot row for partial pivoting
  int pp_row;

  // calls function to perform gaussian pivot with partial pivoting
  for (i = 0; i < m->n - 1; i++) {
    pp_row = find_greatest_in_column(m, i);
    switch_rows(m, B, i, pp_row);
    gaussian_pivot(m, m->n - i);
  }
}

// function that performs the back substitution to solve for our X vector in UX=Y
// takes pointers to structure mat_t (containing matrix A) and vector B as arguments
void UX_equals_Y(mat_t *m, float *B) {

  int i, j;

  for (i = m->n - 1; i >= 0; i--) {
    for (j = m->n - i - 2; j >= 0; j--) {
      B[i] -= m->m[i][m->n - j - 1] * B[m->n - j - 1];
    }
    B[i] /= m->m[i][i];
  }
}

// function that performs the forward substitution to solve for our Y vector in LY=B
// takes pointers to structure mat_t (containing matrix A) and vector B as arguments
void LY_equals_B(mat_t *m, float *B) {

  int i, j;

  for (i = 1; i < m->n; i++) {
    for (j = 0; j < i; j++) {
      B[i] -= m->m[i][j] * B[j];
    }
  }
}

// function to fill sets of test matrices A and test vectors B for algorithm validation
// takes pointers to structure mat_t (containing matrix A) and vector B as arguments, as well as an int to determine which set will be tested
// set 0 would contain a division by 0 if gaussian elimination without partial pivoting were performed
void fill_test_matrices(mat_t *m, float *B, int matrix_set) {

  if (matrix_set == 0) {
    m->m[0][0] = 1;
    m->m[0][1] = -1;
    m->m[0][2] = 3;
    m->m[1][0] = 4;
    m->m[1][1] = -2;
    m->m[1][2] = 1;
    m->m[2][0] = -3;
    m->m[2][1] = -1;
    m->m[2][2] = 4;

    B[0] = 13;
    B[1] = 15;
    B[2] = 8;

    /* 
        Expected X:
                    [ 2
                     -2
                      3]
    */
  }

  if (matrix_set == 1) {
    m->m[0][0] = 1;
    m->m[0][1] = 1;
    m->m[0][2] = -1;
    m->m[1][0] = 1;
    m->m[1][1] = -2;
    m->m[1][2] = 3;
    m->m[2][0] = 2;
    m->m[2][1] = 3;
    m->m[2][2] = 1;

    B[0] = 4;
    B[1] = -6;
    B[2] = 7;

    /* 
        Expected X:
                    [ 1
                      2
                     -1]
    */
  }

  if (matrix_set == 2) {

    m->m[0][0] = -3;
    m->m[0][1] = 2;
    m->m[0][2] = -1;
    m->m[1][0] = 6;
    m->m[1][1] = -6;
    m->m[1][2] = 7;
    m->m[2][0] = 3;
    m->m[2][1] = -4;
    m->m[2][2] = 4;

    B[0] = -1;
    B[1] = -7;
    B[2] = -6;

    /* 
        Expected X:
                    [ 2
                      2
                     -1]
    */
  }

  // this matrix requires partial pivoting during gaussian elimination
  if (matrix_set == 3) {
    m->m[0][0] = 1;
    m->m[0][1] = 1;
    m->m[0][2] = 1;
    m->m[1][0] = 2;
    m->m[1][1] = 2;
    m->m[1][2] = 5;
    m->m[2][0] = 4;
    m->m[2][1] = 6;
    m->m[2][2] = 8;

    B[0] = 1;
    B[1] = 0;
    B[2] = 0;

    /* 
        Expected X:
                    [ 2.333
                     -0.667
                     -0.667]
    */
  }

  // this matrix has rows that are linearly dependent
  if (matrix_set == 4) {
    m->m[0][0] = 1;
    m->m[0][1] = -1;
    m->m[0][2] = 3;
    m->m[1][0] = 2;
    m->m[1][1] = -2;
    m->m[1][2] = 6;
    m->m[2][0] = -3;
    m->m[2][1] = -1;
    m->m[2][2] = 4;

    B[0] = 13;
    B[1] = 15;
    B[2] = 8;
  }

  // this matrix has two identical columns
  if (matrix_set == 5) {
    m->m[0][0] = 1;
    m->m[0][1] = 0;
    m->m[0][2] = 1;
    m->m[1][0] = 2;
    m->m[1][1] = 2;
    m->m[1][2] = 2;
    m->m[2][0] = 1;
    m->m[2][1] = 1;
    m->m[2][2] = 1;

    B[0] = 1;
    B[1] = 1;
    B[2] = 1;
  }

  // this matrix has two identical rows
  if (matrix_set == 6) {
    m->m[0][0] = 1;
    m->m[0][1] = 0;
    m->m[0][2] = 1;
    m->m[1][0] = 2;
    m->m[1][1] = 2;
    m->m[1][2] = 3;
    m->m[2][0] = 1;
    m->m[2][1] = 0;
    m->m[2][2] = 1;

    B[0] = 1;
    B[1] = 1;
    B[2] = 1;
  }

  // this matrix has a column of zeros
  if (matrix_set == 7) {
    m->m[0][0] = 9;
    m->m[0][1] = 0;
    m->m[0][2] = 1;
    m->m[1][0] = 2;
    m->m[1][1] = 0;
    m->m[1][2] = 3;
    m->m[2][0] = 1;
    m->m[2][1] = 0;
    m->m[2][2] = 8;

    B[0] = 1;
    B[1] = 1;
    B[2] = 1;
  }

  // this matrix has a row of zeros
  if (matrix_set == 8) {
    m->m[0][0] = 1;
    m->m[0][1] = 0;
    m->m[0][2] = 7;
    m->m[1][0] = 0;
    m->m[1][1] = 0;
    m->m[1][2] = 0;
    m->m[2][0] = 5;
    m->m[2][1] = 0;
    m->m[2][2] = 1;

    B[0] = 1;
    B[1] = 1;
    B[2] = 1;
  }

  // this matrix has columns that are linearly dependent
  if (matrix_set == 9) {
    m->m[0][0] = 1;
    m->m[0][1] = -1;
    m->m[0][2] = 3;
    m->m[1][0] = 3;
    m->m[1][1] = -8;
    m->m[1][2] = 9;
    m->m[2][0] = -1;
    m->m[2][1] = 8;
    m->m[2][2] = -3;

    B[0] = 13;
    B[1] = 15;
    B[2] = 8;
  }
}

// function to check for duplicate rows or columns in a matrix
// takes a pointer to structure mat_t containing matrix A and ints representing the rows or columns and the axis we are searching in
// returns 1 is duplicate rows or columns are found, 0 if not
int check_for_duplicate(mat_t *m, int row1, int row2, int axis) {

  int j;

  if (axis == 0) {
    for (j = 0; j < m->n; j++) {
      if (m->m[row1][j] != m->m[row2][j]) {
        return 0;
      }
    }
  }

  if (axis == 1) {
    for (j = 0; j < m->n; j++) {
      if (m->m[j][row1] != m->m[j][row2]) {
        return 0;
      }
    }
  }
  return 1;
}

// function to check for zeros in rows or columns in a matrix
// takes a pointer to structure mat_t containing matrix A and ints representing the row or column and the axis we are searching in
// returns 1 if a row or column of zeros is found, 0 if not
int check_for_zeros(mat_t *m, int row, int axis) {

  int i;

  if (axis == 0) {
    for (i = 0; i < m->n; i++) {
      if (m->m[row][i] != 0) {
        return 0;
      }
    }
  }

  if (axis == 1) {
    for (i = 0; i < m->n; i++) {
      if (m->m[i][row] != 0) {
        return 0;
      }
    }
  }

  // printf("\n\nFOUND ZEROES in row/column %d!!\n", row);
  return 1;
}

// function to check for linear dependence in rows or columns in a matrix
// takes a pointer to structure mat_t containing matrix A and ints representing the rows or columns and the axis we are searching in
// returns 1 is linear dependence is found in rows or columns, 0 if not
int check_for_dependence(mat_t *m, int row1, int row2, int axis) {

  int i, zeros_pairs, quotient_flag, zero_in_1, zero_in_2;
  float quotient;

  zeros_pairs = 0;
  quotient_flag = 0;
  zero_in_1 = 0;
  zero_in_2 = 1;

  if (axis == 0) {

    for (i = 0; i < m->n; i++) {
      // if a zero value and non_zero value are found in the same column of two different rows and the inverse is found in another column of those same rows, the rows cannot be dependent
      if (zero_in_1 == 1 && zero_in_2 == 1) {
        return 0;
      }
      // if two rows are composed of all zeros except 1 shared column, the rows are linearly dependent
      if (zeros_pairs == m->n - 1) {
        return 1;
      }
      // a counter for the number of columns in two rows where both are equal to zero
      if (m->m[row1][i] == 0 && m->m[row2][i] == 0) {
        zeros_pairs++;
        continue;
      }
      // activating a flag to indicate we have found a column where the first value is zero and the second is non-zero in the tow rows we are comparing
      if (m->m[row1][i] == 0 && m->m[row2][i] != 0) {
        if (zero_in_1 == 0) {
          zero_in_1 = 1;
        }
        continue;
      }
      // activating a flag to indicate we have found a column where the second value is zero and the first is non-zero in the tow rows we are comparing
      if (m->m[row1][i] != 0 && m->m[row2][i] == 0) {
        if (zero_in_2 == 0) {
          zero_in_2 = 1;
        }
        continue;
      }
      // if we find a column where both rows are non-zero, we obtain the quotient for comparison with other columns and set the flag to search for a second quotient
      if (quotient_flag == 0) {
        quotient = m->m[row1][i] / m->m[row2][i];
        quotient_flag = 1;
        continue;
      }
      // if one quotient has already been obtained, we compare it to the quotient at the current column. If they differ, the rows are not linearly dependent
      if (quotient_flag == 1) {
        if (m->m[row1][i] / m->m[row2][i] != quotient) {
          return 0;
        }
        // if we have made it to the end of the row and no differing quotient have been found, the rows are linearly dependent
        if (i == m->n - 1 && m->m[row1][i] / m->m[row2][i] == quotient) {
          return 1;
        }
      }
    }
  }

  /*  
      This is the same set of calculations as above in this function.
      The difference is that here we are checking columns instead of rows.
  */

  if (axis == 1) {

    for (i = 0; i < m->n; i++) {

      if (zero_in_1 == 1 && zero_in_2 == 1) {
        return 0;
      }

      if (zeros_pairs == m->n - 1) {
        return 1;
      }
      if (m->m[i][row1] == 0 && m->m[i][row2] == 0) {
        zeros_pairs++;
        continue;
      }

      if (m->m[i][row1] == 0 && m->m[i][row2] != 0) {
        if (zero_in_1 == 0) {
          zero_in_1 = 1;
        }
        continue;
      }

      if (m->m[i][row1] != 0 && m->m[i][row2] == 0) {
        if (zero_in_2 == 0) {
          zero_in_2 = 1;
        }
        continue;
      }

      if (quotient_flag == 0) {
        quotient = m->m[i][row1] / m->m[i][row2];
        quotient_flag = 1;
        continue;
      }
      if (quotient_flag == 1) {
        if (m->m[i][row1] / m->m[i][row2] != quotient) {
          return 0;
        }
        if (i == m->n - 1 && m->m[i][row1] / m->m[i][row2] == quotient) {
          return 1;
        }
      }
    }
  }
  return 0;
}

// wrapper function to check for matrix singularity
// takes pointer to structure mat_t containing matrix A as argument
// returns 1 if singularity is found, 0 if not
int check_singularity(mat_t *m) {

  int i, j;

  for (i = 0; i < m->n; i++) {

    if (check_for_zeros(m, i, 0) == 1) {
      printf("\nFound row of zeros at row %d", i);
      return 1;
    }

    if (check_for_zeros(m, i, 1) == 1) {
      printf("\nFound column of zeros at column %d", i);
      return 1;
    }
    for (j = i + 1; j < m->n; j++) {
      if (check_for_duplicate(m, i, j, 0) == 1) {
        printf("\nFound duplicate rows at rows %d and %d", i, j);
        return 1;
      }
      if (check_for_duplicate(m, j, i, 1) == 1) {
        printf("\nFound duplicate columns at columns %d and %d", i, j);
        return 1;
      }
      if (check_for_dependence(m, i, j, 0) == 1) {
        printf("\nThis check\n");
        printf("\nFound linearly dependent rows at rows %d and %d", i, j);
        return 1;
      }
      if (check_for_dependence(m, j, i, 1) == 1) {
        printf("\nFound linearly dependent columns at columns %d and %d", i, j);
        return 1;
      }
    }
  }
  return 0;
}

// wrapper function that runs tests for LU decomposition and solving of systems of linear equations using most aspects of all functions
// takes pointers to structure mat_t (containing matrix A) and vector B
void run_tests(mat_t *m, float *B) {

  int i;

  printf("\nRunning tests now...\n");
  for (i = 0; i < 10; i++) {
    printf("\nTest %d :\n\n", i + 1);
    fill_test_matrices(m, B, i);
    print_mat(m, 0);

    if (check_singularity(m)) {
      printf(", making this a singular matrix with a determinant of 0 and unsuitable for LU decomposition by Gaussian elimination.\n");
      continue;
    }

    print_vec(B, m->n, 0);
    gaussian_elimination(m, B);
    LY_equals_B(m, B);
    UX_equals_Y(m, B);
    print_mat(m, 1);
    print_vec(B, m->n, 1);
    printf("\n");
  }
  printf("\nTests finished successfully.\n\n");
}

// principle wrapper function that regroups the different processes necessary for LU decomposition and solving of systems of linear equations
// takes pointers to structure mat_t (containing matrix A) and vector B
void solve_lu(mat_t *m, float *B) {

  printf("\nSolving for X with LU Decomposition\n\n");
  print_mat(m, 0);
  print_vec(B, m->n, 0);

  if (check_singularity(m)) {
    printf(", making this a singular matrix with a determinant of 0 and unsuitable for LU decomposition by Gaussian elimination.\n");
    printf("\nCome back and try again.\n");
    exit(0);
  }

  gaussian_elimination(m, B);
  LY_equals_B(m, B);
  UX_equals_Y(m, B);
  print_mat(m, 1);
  print_vec(B, m->n, 1);
  printf("\nFinished successfully!\n\n");
}

int main() {

  mat_t m;
  int nb;
  float B[MAXMAT];
  srand(time(NULL));

  // ******************** Comment out this section for a more interactive version ******************************** //

  nb = 3;
  m = create_A(nb);
  create_B(B, nb);
  run_tests(&m, B);

  nb = 10;
  m = create_A(nb);
  create_B(B, nb);
  solve_lu(&m, B);

  // ************************************************************************************************************* //

  

  // ******************** De-comment this section for a more interactive version ******************************** //
  /*
  char choice;

  printf("\nWelcome. Would you first like to run the tests? (Y/N)\n");
  scanf("%c", &choice);

  while (choice != 'y' && choice != 'Y' && choice != 'n' && choice != 'N') {
    printf("\nSorry, I need you to enter a valid answer.\nWould you first like to run the tests? (Y/N)\n");
    scanf("%c", &choice);
  }
  if (choice == 'y' || choice == 'Y') {
    // function call to run tests
    run_tests(&m, B);
    choice = 'n';
  }
  if (choice == 'n' || choice == 'N') {
    printf("\nOkay, so let's get started.\n\nEnter a dimension for the square matrix, any integer from 2-100 :\n");
    scanf("%d", &nb);
  }
  while (nb < 2 || nb > 100) {
    printf("\nSorry, I need you to enter a valid integer input.\nEnter a dimension for the square matrix, any integer from 2-100 :\n");
    scanf("%d", &nb);
  }
  m = create_A(nb);
  create_B(B, nb);
  solve_lu(&m, B);
*/

  return 0;
}