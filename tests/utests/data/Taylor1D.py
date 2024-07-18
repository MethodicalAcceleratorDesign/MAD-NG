#Script to generate the coefficients of the 1D Taylor polynomial (centered in a specific x0) of a generic function up to order n
import sympy as sp

def taylor_series_sinc(x0, n, precision=25):
    # Define the variable and the sinc function
    x = sp.symbols('x')
    sinc = sp.sin(x) / x
    
    # Calculate the Taylor series expansion around x0 up to order n
    taylor_expansion = sp.series(sinc, x, x0, n + 1)
    
    # Extract the coefficients and evaluate them numerically with the specified precision
    coefficients = [taylor_expansion.coeff((x - x0), i).evalf(precision) for i in range(n + 1)]
    
    return coefficients, taylor_expansion

def write_coefficients_to_file(coefficients, filename):
    with open(filename, 'w') as file:
        for coeff in coefficients:
            file.write(f"{coeff}\n")

# Example usage:
x0 = sp.Rational(1,2)  # Expansion around this point
n = 15   # Order of the Taylor series
precision = 25  # Number of decimal places

coefficients, taylor_expansion = taylor_series_sinc(x0, n, precision)
print(f"Taylor series expansion of sinc(x) around x0={x0} up to order {n}:")
print(taylor_expansion)

# Write coefficients to a file
filename = 'your\\path\\coefficients.txt'
write_coefficients_to_file(coefficients, filename)
print(f"Coefficients written to {filename}")
