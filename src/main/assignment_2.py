#Nathan Carney
#COT4500, Assignment 1
import numpy as np
from decimal import Decimal
from scipy.interpolate import CubicSpline

np.set_printoptions(precision=7, suppress=True, linewidth=100)


def neville(x, f, point):
    # x: list of x-coordinates
    # f: list of y-coordinates
    # point: x-coordinate to interpolate
    n = len(x)  # Get the number of data points
    Q = [[0] * n for _ in range(n)]  # Initialize a 2D array of zeros with dimensions n x n
    for i in range(n):
        Q[i][i] = f[i]  # Set the diagonal entries of Q to the corresponding f values
    for k in range(1, n):
        for i in range(n-k):
            j = i + k
            # Compute the polynomial interpolant for x[i:j+1] using Neville's algorithm
            Q[i][j] = ((point - x[j]) * Q[i][j-1] + (x[i] - point) * Q[i+1][j]) / (x[i] - x[j])
    return Q[0][n-1]  # Return the interpolated value at the specified point



def divided_difference_table(x_points, y_points):
    # set up the matrix
    size: int = len(x_points)
    matrix: np.array = np.zeros((size, size), dtype=object)
    
    # fill the matrix
    for index, row in enumerate(matrix):
        row[0] = Decimal(y_points[index])
    
    # populate the matrix (end points are based on matrix size and max operations 
    # we're using)
    for i in range(1, size):
        for j in range(1, i+1):
            # the numerator are the immediate left and diagonal left indices...
            numerator = matrix[i][j-1] - matrix[i-1][j-1]
            # the denominator is the X-SPAN...
            denominator = Decimal(x_points[i]) - Decimal(x_points[i-j])
            operation = numerator / denominator
            matrix[i][j] = operation
    
    # print(matrix)

    # to print out polynomial approximations like the output
    a = np.zeros(3)
    for i in range(1, 4):
        if(i == 1):
            print("[%.15f, " % matrix[i][i], end='')
        elif(i == 3):
            print("%.15f]" % matrix[i][i], end='')
        else:
            print("%.15f, " % matrix[i][i], end='')
    print()

    return matrix


def get_approximate_result(matrix, x_points, value):
    # p0 is always y0 and we use a reoccuring x to avoid having to recalculate x 
    reoccuring_x_span = 1
    reoccuring_px_result = matrix[0][0]
    
    # we only need the diagonals...and that starts at the first row...
    for index in range(1, len(matrix)):
        polynomial_coefficient = matrix[index][index]
        # we use the previous index for x_points....
        reoccuring_x_span *= Decimal((value - x_points[index-1]))
        
        # get a_of_x * the x_span
        mult_operation = polynomial_coefficient * reoccuring_x_span
        # add the reoccuring px result
        reoccuring_px_result += mult_operation
    
    # final result
    return reoccuring_px_result

def apply_div_dif(matrix: np.array):
    size = len(matrix)
    for i in range(2, size):
        for j in range(2, i+2):
            # skip if value is prefilled (we dont want to accidentally recalculate...)
            if j >= len(matrix[i]) or matrix[i][j] != 0:
                continue
            
            # get left cell entry
            left: float = matrix[i][j-1]
            # get diagonal left entry
            diagonal_left: float = matrix[i-1][j-1]
            # order of numerator is SPECIFIC.
            numerator: float = left - diagonal_left
            # denominator is current i's x_val minus the starting i's x_val....
            denominator = matrix[i][0] - matrix[i-(j-1)][0]
            #print("left" , left,"- diagonal left", diagonal_left)
            #print("matrix[i][0]" , matrix[i][0],"- matrix[i][0]", matrix[-][])
            
            # something save into matrix
            operation = numerator / denominator
            matrix[i][j] = operation
    
    return matrix

def hermite_interpolation(x_points, y_points, slope_points):
    # matrix size changes because of "doubling" up info for hermite 
    num_of_points = len(x_points)
    matrix = np.zeros((2*num_of_points, 2*num_of_points))
    # populate x values (make sure to fill every TWO rows)
    for x in range(num_of_points):
        matrix[2*x][0] = x_points[x]
        matrix[2*x+1][0] = x_points[x]
    
    # prepopulate y values (make sure to fill every TWO rows)
    for x in range(num_of_points):
        matrix[2*x][1] = y_points[x]
        matrix[2*x+1][1] = y_points[x]
    
    # # prepopulate with derivates (make sure to fill every TWO rows. starting row CHANGES.)
    for x in range(num_of_points):
        matrix[2*x+1][2] = slopes[x]
    
    filled_matrix = apply_div_dif(matrix)
    print(filled_matrix)

def natural_cubic_spline_matrix(points):
    # Compute the distances hi and the values yi - yi-1 for each interval
    n = len(points) - 1
    h = []
    d = []

    for i in range(n):
        # Compute distance between x-coordinates
        h_i = points[i+1][0] - points[i][0]
        h.append(h_i)
    
        # Compute constants
        d_i = (points[i+1][1] - points[i][1]) / h_i
        d.append(d_i)
    
    # Set up the tridiagonal matrix A
    A = np.zeros((n+1, n+1))
    A[0, 0] = 1
    A[n, n] = 1
    for i in range(1, n):
        A[i, i-1] = h[i-1]
        A[i, i] = 2 * (h[i-1] + h[i])
        A[i, i+1] = h[i]
    
    # Adjust the matrix A for the natural cubic spline
    A[0, 1] = 0
    A[n, n-1] = 0

    b = np.zeros(n+1)
    for i in range(1, n):
        b[i] = 3 * (d[i] - d[i-1])
    b[0] = 0
    b[n] = 0

    # Solve for the vector x
    x = np.linalg.solve(A, b)

    # Return the matrix A
    return A, b, x

if __name__ == "__main__":
    # Neville interpolation
    x = [3.6, 3.8, 3.9]
    f = [1.675, 1.436, 1.318]
    point = 3.7
    result = neville(x, f, point)
    print(result, "\n")

    # Divided difference table
    x_points = [7.2, 7.4, 7.5, 7.6]
    y_points = [23.5492, 25.3913, 26.8224, 27.4589]
    divided_table = divided_difference_table(x_points, y_points)
    print()

    # Approximation
    approximating_x = 7.3
    final_approximation = get_approximate_result(divided_table, x_points, approximating_x)
    print("%.15f" % final_approximation, "\n")

    # Hermite interpolation
    x_points = [3.6, 3.8, 3.9]
    y_points = [1.675, 1.436, 1.318]
    slopes = [-1.195, -1.188, -1.182]
    hermite_interpolation(x_points, y_points, slopes)
    print()

    # Natural cubic spline matrix
    points = [(2,3), (5,5), (8,7), (10,9)]
    A, b, x = natural_cubic_spline_matrix(points)
    print(A, "\n")
    print(b, "\n")
    print(x, "\n")

   
