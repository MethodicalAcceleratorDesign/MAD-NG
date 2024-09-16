import sympy as sp

def compute_and_evaluate_derivatives(x0_points, n, fun_expr, precision=25):

    all_derivatives = []

    for x0 in x0_points:
        derivatives_at_x0 = []
        for i in range(n + 1):
            # Compute the i-th derivative of the function
            derivative = sp.diff(fun_expr, x, i)
            
            # Divide by factorial(i)
            derivative /= sp.factorial(i)
            
            # Evaluate the derivative at x0
            derivative_at_x0 = derivative.subs(x, x0).evalf(precision + 5)  # Extra precision
            derivatives_at_x0.append(derivative_at_x0)
        
        all_derivatives.append(derivatives_at_x0)
    
    return all_derivatives

def format_x0_as_pi_notation(x0):
    # Check if x0 is a multiple of pi
    if sp.simplify(x0 / sp.pi).is_rational:
        coeff = sp.simplify(x0 / sp.pi)
        return f"pi*{coeff}".replace('*1/', '/') if coeff != 1 else "pi"
    else:
        # For complex numbers, handle separately
        real_part = sp.re(x0)
        imag_part = sp.im(x0)
        
        real_str = (f"pi*{sp.simplify(real_part / sp.pi)}" if sp.simplify(real_part / sp.pi).is_rational else str(real_part.evalf())).replace('*1/', '/')
        imag_str = (f"pi*{sp.simplify(imag_part / sp.pi)}" if sp.simplify(imag_part / sp.pi).is_rational else str(imag_part.evalf())).replace('*1/', '/')
        
        if real_part == 0:
            return f"{imag_str}*1i"
        elif imag_part == 0:
            return real_str
        else:
            return f"{real_str} + {imag_str}*1i" if imag_part > 0 else f"{real_str} - {imag_str[1:]}*1i"

def write_coefficients_to_file(all_derivatives, x0_points, filename, function_name):
    with open(filename, 'w') as file:
        # Write the header with the function name
        file.write(f"M.fun.{function_name} = {{\n")
        
        for i, derivatives in enumerate(all_derivatives):
            x0 = x0_points[i]
            
            # Format x0 in pi notation if possible
            x0_str = format_x0_as_pi_notation(x0)
            file.write(f"{{ x0 = {x0_str},\n")
            
            for j, derivative in enumerate(derivatives):
                if j % 2 == 0 and j != 0:
                    file.write("\n")
                
                # Check if the derivative is real or complex
                if derivative.is_real:
                    # Write the real derivative in scientific notation
                    derivative_str = f"{float(derivative.evalf()):.25e}".replace('e+', 'e+0').replace('e-', 'e-0')
                else:
                    # Format the complex derivative
                    real_part = sp.re(derivative).evalf()
                    imag_part = sp.im(derivative).evalf()

                    real_str = f"{float(real_part):.25e}".replace('e+', 'e+0').replace('e-', 'e-0')
                    imag_str = f"{float(imag_part):.25e}".replace('e+', 'e+0').replace('e-', 'e-0')

                    if real_part != 0:
                        derivative_str = f"{real_str} + {imag_str}i" if imag_part >= 0 else f"{real_str} - {imag_str[1:]}i"
                    else:
                        derivative_str = f"{imag_str}i" if imag_part >= 0 else f"-{imag_str[1:]}i"
                
                file.write(f"  {derivative_str},")
            
            file.write("\n},\n\n")
        
        file.write("}\n")

# Example usage:
if __name__ == "__main__":
    # Set of x0 points
    x0_points = [
        sp.Rational(1, 1000000),
        sp.Rational(1000001, 1000000),
        sp.Rational(1, 1000000) * sp.I,
        sp.Rational(1, 10) * sp.I,
        sp.Rational(1, 1000000) + 1 * sp.I,
        sp.Rational(1000001, 1000000) * sp.I,
        sp.pi / 3 * sp.I,
        2 * sp.I,
        1 + 1 * sp.I,
        10 + 10 * sp.I,
    ]

    # Order of derivatives
    n = 15

    # Number of decimal places
    precision = 25

    x = sp.symbols('x')

    #put here the function you want to expand
    fun_expr = sp.atanh(x)

    # Compute and evaluate the derivatives at the specified points
    all_derivatives = compute_and_evaluate_derivatives(x0_points, n, fun_expr, precision)

    # Output filename
    output_filename = 'your\\path\\function_derivatives_at_points.txt'

    # Specify the name of the function (cannot be automatized because some functions we use either they have different names or they do not exist natively in sympy)
    function_name = "sin"  

    # Write all derivatives to the file
    write_coefficients_to_file(all_derivatives, x0_points, output_filename, function_name)
    print(f"All derivatives written to {output_filename}")

    print("All files generated successfully.")
