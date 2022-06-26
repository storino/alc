from tkinter import *
import json
import copy

def dot(v1, v2):
    return sum([i*j for (i, j) in zip(v1, v2)])

def simetrica(a,n):
    for i in range(n):
        for j in range(n):
            if a[i][j] != a[j][i]:
                return False
    return True            
            
def norma(v):
    return (sum([v[i]**2 for i in range(len(v))]))**0.5

def transpose(a, n):
    for i in range(1,n):
        for j in range(i):
            a[i][j],a[j][i] = a[j][i],a[i][j]

def diag_dominante(a, n):
    linha,coluna = 1,1
    for i in range(n):
        if sum(a[i]) > 2*a[i][i]:
            linha = 0
    transpose(a,n)
    for i in range(n):
        if sum(a[i]) > 2*a[i][i]:
            coluna = 0
    return linha or coluna

def jacobi(a, b, n, tol=10**(-5)):
    x0 = [1.0 for i in range(n)]
    x1 = [0.0 for i in range(n)]
    r = 1
    residuos = []
    while r > tol:
        for i in range(n):
            soma = sum([a[i][j]*x0[j] if i != j else 0 for j in range(n)])
            x1[i] = (b[i]-soma)/a[i][i]
        r = round(norma([i-j for (i, j) in zip(x1, x0)])/norma(x1), 6)
        residuos.append(r)
        x0 = x1[:]
    x1 = [round(i,6) for i in x1]
    return x1, residuos, len(residuos)

def seidel(a, b, n, tol=10**(-5)):
    if n != len(a[0]):
        print("Matriz não é quadrada")
        return
    if (not diag_dominante(a,n)):
        print("Não é diagonal dominante")
        return
    x0 = [1.0 for i in range(n)]
    x1 = [0.0 for i in range(n)]
    r = 1
    residuos = []
    while r > tol:
        for i in range(n):
            soma1 = sum([a[i][j]*x1[j] for j in range(i)])
            soma2 = sum([a[i][j]*x0[j] for j in range(i+1,n)])
            x1[i] = (b[i]-soma1-soma2)/a[i][i]
        r = round(norma([i-j for (i, j) in zip(x1, x0)])/norma(x1), 6)
        residuos.append(r)
        x0 = x1[:]
    x1 = [round(i,6) for i in x1]
    return x1, residuos, len(residuos)

def LUdecomp(a, n):
    for k in range(0,n-1):
        for i in range(k+1,n):
            if a[i][k] != 0.0:
                m = a[i][k] / a[k][k]
                for j in range(k+1,n):
                    a[i][j] = a[i][j] - m*a[k][j]
                a[i][k] = m
    return a

def cholesky(a, n):
    g = [[0.0]*n for i in range (n)]
    pd = 1
    if a[0][0] <= 0: pd = 0
    g[0][0] = (a[0][0])**0.5
    for i in range(1,n):
        g[i][0]=a[i][0]/g[0][0]
    for j in range(1,n):
        for i in range(j,n):
            soma = 0
            if j==i:
                soma = sum(g[i][k]**2 for k in range(i))
                if a[i][i] - soma <= 0: pd = 0
                g[i][i] = (a[i][i] - soma)**0.5
            else:
                soma = sum(g[i][k]*g[j][k] for k in range(j))
                g[i][j] = (a[i][j] - soma)/g[j][j]
    return g, pd

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

def auxiliar(a, index):
    secondary = copy.deepcopy(a)

    for row in range(len(secondary)):
        secondary[row] = secondary[row][1:]

    return secondary[:index] + secondary[index+1:]

def solver(a, b, ICOD, n, tol=10**(-5)):
    if ICOD == 1:
        a = LUdecomp(a,n)
        for k in range(1,n):
            b[k] = b[k] - dot(a[k][0:k],b[0:k])
        b[n-1] = b[n-1]/a[n-1][n-1]
        for k in range(n-2,-1,-1):
            b[k] = (b[k] - dot(a[k][k+1:n],b[k+1:n]))/a[k][k]
        return b
    elif ICOD == 2:
        if (not simetrica(a,n)):
            return "A matriz não é simétrica"
        a, pd = cholesky(a,n)
        if pd == 0: return "A matriz não é positiva definida"
        y = [0.0 for i in range(n)]
        for i in range(n):
            soma = sum([a[i][j]*y[j] for j in range(i)])
            y[i] = (b[i]-soma)/a[i][i]
        x = [0.0 for i in range(n)]
        transpose(a,n)
        x[n-1]=round(y[n-1]/a[n-1][n-1],5)
        for i in range(n-2,-1,-1):
            soma = sum([a[i][j]*x[j] for j in range(i+1,n)])
            x[i] = round((y[i]-soma)/a[i][i],5)
        return x
    elif ICOD == 3:
        if (not diag_dominante(a,n)):
            return "A matriz não é diagonal dominante"
        return jacobi(a, b , n, tol)
    elif ICOD == 4:
        return seidel(a, b, n, tol)
    else:
        return "ICOD inválido"

root = Tk()
root.title("task1")
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
i_A.insert(0,"[[3,-1,-1],[-1,3,-1],[-1,-1,3]]")
i_A.grid(row=3, column=1)

i_b = Entry(root, width=100)
Label(root, text="b : ").grid(row=4, column=0)
i_b.insert(0, "[1,2,1]")
i_b.grid(row=4, column=1)

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
    b = json.loads(i_b.get())
    Det = determinant(A)
    if ICOD == 3 or ICOD == 4:
        TOLm = json.loads(i_TOLm.get())
        result = solver(A, b, ICOD, n, TOLm)
        if type(result) == str: r = result
        else:
            r = f"X = {result[0]}\n Resíduos = {result[1]}\n Nº de Iterações = {result[2]}"
    elif ICOD == 1 or ICOD == 2:
        result = solver(A, b, ICOD, n)
        if type(result) == str: r = result
        else:
            r = f"X = {result}"
    if IDET > 0:
        r += f"\n Det = {Det}"
    myLabel = Label(root, text = r)
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