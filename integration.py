import math
import cmath
import numpy as np
from scipy.integrate import quad

#########################################################################################################
########################################The calculation of integration################################### 
#########################################################################################################
# a)  calculate the following integration:
    
# Define the function f(x)
def f(x):
    return math.cos(x) + 6j * math.sin(2*x)
'''
def f_times_exp(x, w):
    return  cmath.exp(1j * x * w)
'''
# Define the trapezoidal rule for numerical integration
def trapezoidal_rule(a, b, n, w ,f):
    h = (b - a) / n
    integral = 0
    for i in range(n+1):
        x = a + i * h
        if i == 0 or i == n:
            integral += f(x) * cmath.exp(1j * x * w) / 2
        else:
            integral += f(x) * cmath.exp(1j * x * w)
    integral *= h
    return integral
# Input from the user
a = float(input("Enter the lower limit (a): "))
b = float(input("Enter the upper limit (b): "))
w = float(input("Enter the value of w: "))
n = int(input("Enter the number of subintervals: "))
# Compute the integral using the trapezoidal rule
integral_value = trapezoidal_rule(a, b, n, w,f)
print("The value of the integral is:", integral_value)

#########################################################################################################
# b) show the results in both Cartesian and polar forms:
    
# Function to convert Cartesian to polar form
def cartesian_to_polar(z):
    r = abs(z)
    theta = np.angle(z)  # Use numpy's angle function for correct phase calculation
    return r, theta
# Convert results to Cartesian and polar forms
cartesian_result = integral_value.real, integral_value.imag
polar_result = cartesian_to_polar(integral_value)
print("Integration Result (Cartesian Form):", cartesian_result)
print("Integration Result (Polar Form):", polar_result)

#########################################################################################################
# c) Compare you results with respect to the analytical form:
    
# Define the antiderivative function F(x)
def antiderivative_F(x):
    return cmath.sin(x) - 3j * cmath.cos(2*x)
# Analytical form of the integral
def analytical_integral(a, b, w):
    F_b = antiderivative_F(b)
    F_a = antiderivative_F(a)
    return (F_b - F_a) * cmath.exp(1j * b * w)
# Compute the analytical integral
analytical_result = analytical_integral(a, b, w)
print("The analytical form of the integral is:", analytical_result)

#########################################################################################################
# Using built-in function

def integrand(x, w):
    return f(x) * np.exp(1j * w * x)
result_real, error_real = quad(lambda x: np.real(integrand(x, w)), a, b)
result_imag, error_imag = quad(lambda x: np.imag(integrand(x, w)), a, b)

print("Using Built-in Function:")
print("Real part:", result_real)
print("Imaginary part:", result_imag)

#########################################################################################################
# d) apply the same code to functions different from the above f(x) to contain two different combinations:

# Define the first function f(x)
def f1(x):
    return math.cos(x) + math.log(x**2)
integral_value_f1 = trapezoidal_rule(a, b, n, w, f1)
print("The value of the integral for f1(x) is:", integral_value_f1)

# Define the second function f(x)
def f2(x):
    return math.tan(x) + math.exp(x)
integral_value_f2 = trapezoidal_rule(a, b, n, w, f2)
print("The value of the integral for f2(x) is:", integral_value_f2)
#########################################################################################################

