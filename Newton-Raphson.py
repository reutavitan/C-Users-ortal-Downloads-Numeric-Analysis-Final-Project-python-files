import sympy as sp
from sympy.utilities.lambdify import lambdify
from datetime import datetime


# Defining Function
def f(x):
    return (x * sp.exp(-x**2 +5*x)) * (2 * x ** 2 - 3 * x - 5)

# Defining derivative of function
def g(z):
    x = sp.symbols('x')
    f = (sp.sin((2 * x ** 3) + (5 * x ** 2) - 6)) / (2 * sp.exp((-2) * x))
    f_prime = f.diff(x)  # Derivation of  f by x
    f_prime = lambdify(x, f_prime)
    return f_prime(z)


# Implementing Newton Raphson Method
def newtonRaphson(x0, e):
    print('\n\n*** NEWTON RAPHSON METHOD IMPLEMENTATION ***')
    step = 1
    condition = True
    while condition:
        if g(x0) == 0.0:
            print('Divide by zero error!')
            break

        x1 = x0 - f(x0) / g(x0)
        print('Iteration-%d, x1 = %0.6f and f(x1) = %0.6f' % (step, x1, f(x1))+"00000"+d+h+m)
        x0 = x1
        step = step + 1

        condition = abs(f(x1)) > e

    print('\nRequired root is: %0.8f' % x1+"00000"+d+h+m)

local_dt = datetime.now()
d = str(local_dt.day)
h = str(local_dt.hour)
m = str(local_dt.minute)
# Input Section
# Function (x * sp.exp(-x**2 +5*x)) * (2 * x ** 2 - 3 * x - 5)
x1 = float(input('Enter start Range: '))
x2 = float(input('Enter end Range: '))
x0 = 0.1
# Converting x0
x0 = float(x0)
e = float(0.00001)
# Starting Newton Raphson Method
newtonRaphson(x0, e)
import math
from datetime import datetime
import sympy as sp
from sympy.utilities.lambdify import lambdify
x = sp.symbols('x')
my_f = (x * sp.exp(-x**2 +5*x)) * (2 * x ** 2 - 3 * x - 5)
def SecantMethodInRangeIterations(f, check_range, epsilon=0.0001):
    roots = []
    iterCounter = 0
    for i in check_range:
        startPoint = round(i, 2)
        endPoint = round(i + 0.1, 2)
        print("Checked range:", startPoint, "-", endPoint)
        # Send to the Secant Method with 2 guesses
        local_root = SecantMethod(f, startPoint, endPoint, epsilon, iterCounter)
        # If the root has been found in previous iterations
        if round(local_root, 6) in roots:
            print("Already found that root.")
        # If the root is out of range tested
        elif not (startPoint <= local_root <= endPoint):
            print("root out of range.")
        elif local_root is not None:
            roots += [round(local_root, 6)]
    return roots
def SecantMethod(func, firstGuess, secondGuess, epsilon, iterCounter):
    if iterCounter > 100:
        return

    if abs(secondGuess - firstGuess) < epsilon:  # Stop condition
        print("after ", iterCounter, "iterations The root found is: ", round(secondGuess, 6) , "00000" + d + h + m)
        return round(secondGuess, 6)  # Returns the root found

    next_guess = (firstGuess * func(secondGuess) - secondGuess * func(firstGuess)) / (
                func(secondGuess) - func(firstGuess))
    print("iteration no.", iterCounter, "\tXi = ", firstGuess, " \tXi+1 = ", secondGuess,
          "\tf(Xi) = ", func(firstGuess))
    # Recursive call with the following guess
    return SecantMethod(func, secondGuess, next_guess, epsilon, iterCounter + 1)
local_dt = datetime.now()
d = str(local_dt.day)
h = str(local_dt.hour)
m = str(local_dt.minute)
def frange(start, end=None, inc=None):
    "Function for a range with incomplete numbers"
    if end == None:
        end = start + 0.0
        start = 0.0

    if inc == None:
        inc = 1.0

    L = []
    while 1:
        next = start + len(L) * inc
        if inc > 0 and next >= end:
            break
        elif inc < 0 and next <= end:
            break
        L.append(next)

    return L
checkRange = frange(0, 3, 0.1)
def func(val):
    return lambdify(x, my_f)(val)



epsilon = 0.0001
# Function (x * sp.exp(-x**2 +5*x)) * (2 * x ** 2 - 3 * x - 5)
checkRange = frange(-1, 1.6, 0.1)
print("\n*** Secant Method***")
SecantMethodInRangeIterations(func, checkRange, 0.0001)
import math
import sympy as sp
from sympy.utilities.lambdify import lambdify

x = sp.symbols('x')
my_f = (x * sp.exp(-x**2 +5*x)) * (2 * x ** 2 - 3 * x - 5)


def printFinalResult(result):
    """
    Function for printing final results according to the requested format
    :param result: The final results (list or number)
    :return: print the result
    """
    from datetime import datetime
    local_dt = datetime.now()
    d = str(local_dt.day)
    h = str(local_dt.hour)
    m = str(local_dt.minute)
    if isinstance(result, list):
        for i in range(len(result)):
            print(str(result[i]) + "00000" + d + h + m)
        return
    return str(result) + "00000" + d + h + m


def SimpsonRule(func, n, a, b):
    if n % 2 != 0:
        return 0, False
    h = (b - a) / n
    print("h = ", h)
    str_even = ""
    str_odd = ""
    k2 = b
    # print("Error evaluation En = ", round(SimpsonError(my_f, b, a, h), 6))
    integral = func(a) + func(b)
    # Calculation of a polynomial lagranz for the sections
    for i in range(n):
        k = a + i * h  # new a
        if i != n - 1:  # new b
            k2 = a + (i + 1) * h
        if i % 2 == 0:  # even places
            integral += 2 * func(k)
            str_even = "2 * " + str(func(k))
        else:  # odd places
            integral += 4 * func(k)
            str_odd = "4 * " + str(func(k))
        print("h/3 ( ", str(func(k)), " + ", str_odd, " + ", str_even, " + ", str(func(k2)), " )")
    integral *= (h / 3)
    return integral, True


epsilon = 0.0001
n = 20


def func(val):
    return lambdify(x, my_f)(val)


print("\n***  Simpson’s  ***")
res = SimpsonRule(func, n, 0.5, 1)
if res[1]:
    print("I = ", printFinalResult(round(res[0], 6)))
else:
    print("n must be even !")
