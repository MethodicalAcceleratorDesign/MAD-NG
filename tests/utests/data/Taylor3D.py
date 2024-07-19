import sympy as sp
import numpy as np

def sinc_function(x):
    return sp.sinh(x) / x

def taylor_series_sinc(x0, n):
    # Define the variable and the sinc function
    x = sp.symbols('x')
    sinc = sp.sinh(x) / x
    
    # Calculate the Taylor series expansion around x0 up to order n
    taylor_expansion = sp.series(sinc, x, x0, n + 1)
    
    # Extract the coefficients symbolically
    coefficients = [taylor_expansion.coeff((x - x0), i) for i in range(n + 1)]
    
    # Replace the first coefficient with the function value at x0
    if x0 == 0: 
        coefficients[0] = sp.sympify(1)
    else:
        coefficients[0] = sinc_function(x0)
    
    return coefficients

def read_coefficients(filename):
    """Reads the coefficients and indices from the file and returns them as a list of tuples."""
    coefficients = []
    with open(filename, 'r') as file:
        for line in file:
            parts = line.split()
            if len(parts) == 2:
                coefficient = int(parts[0])
                index = int(parts[1])
                coefficients.append((coefficient, index))
    return coefficients

def multiply_coefficients_with_vector(coefficients, vector):
    """Multiplies each coefficient with the value at the specified index in the vector symbolically."""
    new_vector = []
    for coeff, idx in coefficients:
        new_value = coeff * vector[idx]
        new_vector.append(new_value)
    return new_vector

def write_coefficients_to_file(coefficients_list, x0_points, filename, columns=4):
    with open(filename, 'w') as file:
        file.write("M.fun.asinhc = {\n")
        
        for i, coefficients in enumerate(coefficients_list):
            x0 = x0_points[i]
            file.write(f"{{ x0 = {float(x0)},\n")
            
            for j, coeff in enumerate(coefficients):
                if j % columns == 0 and j != 0:
                    file.write("\n")
                file.write(f"  {float(coeff):.20e},")
                
            file.write("\n},\n\n")
        
        file.write("}\n")

# Set of x0 points
x0_points = [
    0,
    sp.Rational(1,1000),
    sp.Rational(1,10),
    1,
    10,
    sp.pi/5,
    sp.pi/3
]

# Order of the Taylor series
n = 15

# Generate coefficients for each x0
coefficients_list = []
for x0 in x0_points:
    coefficients = taylor_series_sinc(x0, n)
    coefficients_list.append(coefficients)

# Read multinomial coefficients from file
multinomial_filename = 'C:\\Users\\anton\\Downloads\\coefficients.txt'
multinomial_coeffs = read_coefficients(multinomial_filename)

# Apply multinomial coefficients to each set of Taylor series coefficients and evaluate using SymPy
adjusted_coefficients_list = []
for coeffs in coefficients_list:
    adjusted_coeffs = multiply_coefficients_with_vector(multinomial_coeffs, coeffs)
    # Evaluate each coefficient using SymPy
    evaluated_coeffs = [sp.N(coeff, 25) for coeff in adjusted_coeffs]
    adjusted_coefficients_list.append(evaluated_coeffs)

# Output filename
output_filename = 'C:\\Users\\anton\\Downloads\\formatted_coefficients.txt'

# Write all adjusted coefficients to the same file
write_coefficients_to_file(adjusted_coefficients_list, x0_points, output_filename, columns=4)
print(f"All formatted coefficients written to {output_filename}")

print("All files generated successfully.")
