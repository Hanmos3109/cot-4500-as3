import numpy

numpy.set_printoptions(precision=7, suppress=True, linewidth=100)

def func(t: float, w: float):
    return t - pow(w,2)


def eulerMet(r, n, init):
    h = (r[1] - r[0])/n # .2
    i = 0
    w_0 = init[1]
    w_current = w_0
    for x in range(10):
        t_i = r[0] + h*i #0 + .2*i  (0, .2, .4, ...)
        w_next = w_current + h*func(t_i, w_current)
        
        
        w_current = w_next
        i = i+1
    print(w_current)
    print("")

def rungeKutta(r, n, init):
    h = (r[1] - r[0])/n # .2
    i = 0
    w_0 = init[1]
    w_current = w_0
    for x in range(10):
        t_i = r[0] + h*i #0 + .2*i  (0, .2, .4, ...)
        t_i1 = r[0] + h*(i+1)
        k1 = h*func(t_i, w_current)
        k2 = h*func(t_i+(h/2), w_current + .5*k1)
        k3 = h*func(t_i+(h/2), w_current + .5*k2)
        k4 = h*func(t_i1, w_current+k3)

        w_next = w_current + (1/6)*(k1+2*k2+2*k3+k4)
        
        
        w_current = w_next
        i = i+1
    print(w_current) 
    print("")

def gauss(matrix):
    n = 3
    c = 1
    e1 = matrix[0]
    
    for i in range(1, 3):
        matrix[c] = (matrix[c] - ((matrix[c][0])/(matrix[0][0]))*e1)
        
        c = c + 1
    
    c = 2
    for i in range(2, 3):
        matrix[c] = (matrix[c] - ((matrix[c][1])/(matrix[1][1]))*matrix[c-1])      
        
        c = c + 1
    

    xArray = numpy.zeros(3)
    
    t = n-2
    xArray[n-1] = matrix[(t+1)][t+2]/matrix[t+1][t+1]
    xArray[n-2] = (matrix[t][t+2]-matrix[t][t+1]*xArray[t+1])/(matrix[t][t])
    xArray[n-3] = (matrix[t-1][t+2]-matrix[t-1][t+1]*xArray[t+1]-matrix[t-1][t]*xArray[t])/(matrix[t-1][t-1])
    
    print(xArray)

    print("")



def LUMatrix(matrix):
    #Determinant
    print(numpy.linalg.det(matrix))
    print("")
    n = 4
    c1 = 0
    c2 = 0
    lmatrix = numpy.zeros((n, n))
    for i in range(0, n):
        #Setting up Identity Matrix
        c2 = 0
        for j in range(0, n):
            if c1 == c2:
                lmatrix[c1][c2] = 1
            c2 = c2 + 1
        c1 = c1 + 1


    #
    #UPPER TRIANGLE
    #
    
    r = 0
    c = 1
    count = 1
    for i in range(1, n): #R
        c = 1
        for j in range(1, n): #C
            if c >= count:
                lmatrix[c][r] = matrix[c][r]/matrix[r][r]
                matrix[c] = (matrix[c] - ((matrix[c][r])/(matrix[r][r]))*matrix[r])
                
                     
            c = c + 1 
        r = r + 1
        count = count + 1
    print(lmatrix) #LOWER
    print("")
    print(matrix) #UPPER
    print("")

def isDD(matrix):
    c = 0
    r = 0
    n = 4 #(0 -> 4)
    flag = "True"
    diagArray = numpy.zeros(n+1)
    for i in range(0,n+1): #Rows
        c = 0
        sum = 0
        for j in range(0,n+1): #Columns
            if c == r:
                diagArray[c] = abs(matrix[c][r])
            else:
                sum = sum + abs(matrix[c][r])
            c = c + 1 
        if sum > diagArray[r]:
            flag = "False"
        r = r + 1

    if flag == "False":
        print("False")
    else:
        print("True")
    print("")

def isPD(matrix):
    symFlag = "True"
    posFlag = "True"
    n = 3
    #Check is symetrical
    c = 0
    r = 0
    for i in range(0, n):
        c = 0
        for j in range(0, n):
            if c > r:
                if matrix[c][r] != matrix[r][c]:
                    symFlag = "False"
                    posFlag = "False"
            c = c + 1
        r = r + 1
    if symFlag == "True":
        x = numpy.linalg.eigh(matrix)
        c1 = 0
        for i in range(0, n):
            if x[0][c1] < 0:
                posFlag = "False"
            c1 = c1 + 1

    
    print(posFlag)

    





# Main body of program
# Problem 1
# Euler's Method
funcRange = [0,2]
numIters = 10
initPoint = [0,1]
eulerMet(funcRange, numIters, initPoint)

# Problem 2
# Runge-Kutta
funcRange = [0,2]
numIters = 10
initPoint = [0,1]
rungeKutta(funcRange, numIters, initPoint)

# Problem 3 
# Gaussian Elimination with backwards substitution

aMatrix = numpy.array([[ 2, -1,  1,  6],
                       [ 1,  3,  1,  0],
                       [-1,  5,  4, -3]], dtype=numpy.double)
gauss(aMatrix)

# Problem 4
# LU Factorization
bMatrix = numpy.array([[ 1,  1,  0,  3],
                       [ 2,  1, -1,  1],
                       [ 3, -1, -1,  2],
                       [-1,  2,  3, -1]])

LUMatrix(bMatrix)

# Problem 5
# Diagonally Dominant?
cMatrix = numpy.array([[ 9,  0,  5,  2,  1],
                       [ 3,  9,  1,  2,  1],
                       [ 0,  1,  7,  2,  3],
                       [ 4,  2,  3, 12,  2],
                       [ 3,  2,  4,  0,  8]])
isDD(cMatrix)

# Problem 6
# Positive Definite?
dMatrix = numpy.array([[ 2,  2,  1],
                       [ 2,  3,  0],
                       [ 1,  0,  2]])
isPD(dMatrix)


