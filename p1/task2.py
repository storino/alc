import math
from tkinter import *
import json
import copy

def auxiliar(a, index):
    secondary = copy.deepcopy(a)

    for row in range(len(secondary)):
        secondary[row] = secondary[row][1:]

    return secondary[:index] + secondary[index+1:]

def simetrica(a,n):
    for i in range(n):
        for j in range(n):
            if a[i][j] != a[j][i]:
                return False
    return True

def determinant(a):
    number_of_rows = len(a)
    number_of_columns = len(a[0])
    result = 0
    if(number_of_rows != number_of_columns):
        return "Não é possível calcular a determinante de uma matriz não quadrada"

    if(number_of_rows == 1):
        return a[0][0]

    for k in range(number_of_rows):
        result += a[k][0]*((-1)**k) * \
            determinant(auxiliar(a, k))

    return result

def transpose(a, n):
    b = copy.copy(a)
    for i in range(1,n):
        for j in range(i):
            b[i][j],b[j][i] = b[j][i],b[i][j]
    return b

def mult_matrix_v(a, v):
    n = len(a)
    m = range(len(v))
    result = [0.0 for i in range(n)]
    for j in range(n):
        result[j] = sum([a[j][i]*v[i] for i in m])
    return result

def mult_matrix(a, b):
    n = len(a)
    result = [[0.0 for i in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(len(b[0])):
            for k in range(len(b)):
                result[i][j] += a[i][k]*b[k][j]
    return result

def potencia(a, n, tol=10**(-5)):
    x = [1.0 for i in range(n)]
    lam1 = 1.0
    residuos = []
    r = 1
    while r > tol:
        lam0 = lam1
        y = mult_matrix_v(a,x)
        lam1 = y[0]
        x = [i/lam1 for i in y]
        r = abs(lam1-lam0)/abs(lam1)
        residuos.append(round(r,5))
    lam1 = round(lam1,5)
    x = [round(i,5) for i in x]
    return lam1, x, residuos, len(residuos)

def maior_elemento(a, n):
    maior = 0
    for i in range(n):
        for j in range(n):
            if (i != j) and ( float(abs(a[i][j])) > maior ):
                maior = float(abs(a[i][j]))
                indices = (i , j)
    return indices

def calculate_phi(a, indices):
    r, s = indices
    if a[r][r] == a[s][s]:
        return math.pi/4
    else:
        return math.atan(2*a[r][s]/(a[r][r]-a[s][s])) / 2

def calculate_p(a, n, indices):
    r,s = indices
    phi = 0
    p = [[float(i == j) for j in range(n)] for i in range(n)]
    phi = calculate_phi(a, indices)
    p[r][r] = math.cos(phi)
    p[s][s] = math.cos(phi)
    p[r][s] = -math.sin(phi)
    p[s][r] = math.sin(phi)
    return p

def jacobi(a, n , tol=10**(-5)):
    r,s = maior_elemento(a, n)
    x =[[float(i == j) for j in range(n)] for i in range(n)]
    count = 0
    while float(abs(a[r][s])) > tol:
        p = calculate_p(a, n, (r,s) )
        temp = mult_matrix(a,p)
        p_t = transpose(p, n)
        a = mult_matrix(p_t,temp)
        x = mult_matrix(x,p)
        r,s = maior_elemento(a, n)
        count += 1
    autovalores = []
    autovetores = []
    for i in range(n):
        autovalores.append(round(a[i][i],5))
        autovetores.append([round(x[i][j],5) for j in range(n)])
    return autovalores, autovetores, count

def solver(a,n,ICOD,TOLm):
    if ICOD == 1:
        x = potencia(a, n, TOLm)
        return f"Maior autovalor = {x[0]} \n Autovetor associado = {x[1]} \n Resíduos = {x[2]} \n Iterações = {x[3]}"
    if ICOD == 2:
        x = jacobi(a, n , TOLm)
        return f"Autovalores = {x[0]} \n Autovetores = {x[1]} \n Iterações = {x[2]}"

root = Tk()
root.title("task2")
i_n = Entry(root, width=100)
Label(root, text="N : ").grid(row=0, column=0)
i_n.insert(0, "3")
i_n.grid(row=0, column=1)

i_ICOD = Entry(root, width=100)
Label(root, text="ICOD : ").grid(row=1, column=0)
i_ICOD.insert(0, "1")
i_ICOD.grid(row=1, column=1)

i_IDET = Entry(root, width=100)
Label(root, text="IDET : ").grid(row=2, column=0)
i_IDET.insert(0,"1")
i_IDET.grid(row=2, column=1)

i_A = Entry(root, width=100)
Label(root, text="A : ").grid(row=3, column=0)
i_A.insert(0,"[[1.0,0.2,0.0],[0.2,1.0,0.5],[0.0,0.5,1.0]]")
i_A.grid(row=3, column=1)

i_TOLm = Entry(root, width=100)
Label(root, text="TOLm : ").grid(row=5, column=0)
i_TOLm.insert(0, "0.001")
i_TOLm.grid(row=5, column=1)

def myClick():
    global myLabel
    n = json.loads(i_n.get())
    ICOD = json.loads(i_ICOD.get())
    IDET = json.loads(i_IDET.get())
    A = json.loads(i_A.get())
    TOLm = json.loads(i_TOLm.get())
    Det = determinant(A)
    result = solver(A, n, ICOD, TOLm)
    if ICOD ==2:
        if (not simetrica(A,n)):
            result = "A matriz não é simétrica"
    if IDET > 0:
        result += f"\n Det = {Det}"
    myLabel = Label(root, text = result)
    myLabel.grid(row = 7, column = 1)
    myButton['state'] = DISABLED

def Deletar():
    myLabel.destroy()
    myButton['state'] = NORMAL

myButton = Button(root, text = "Calcular", command = myClick)
myButton.grid(row=6, column=1)

myButton2 = Button(root, text = "Resetar", command = Deletar)
myButton2.grid(row=8, column=1)

root.mainloop()