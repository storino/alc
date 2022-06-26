import numpy as np
from tkinter import Label, Tk, Button, NORMAL, DISABLED, Entry
import json

def function(c, theta):
    N = len(c)
    f = [0 for i in range(N)]
    f[0] = 2*c[1]**2 + c[0]**2 + 6*c[2]**2 - 1.0
    f[1] = 8*c[1]**3 + 6*c[1]*c[0]**2 + 36*c[0]*c[1]*c[2] + \
           108*c[1]*c[2]**2 - theta[0]
    f[2] = 60*c[1]**4 + 60*c[1]**2*c[0]**2 + 576*c[1]**2*c[0]*c[2] + \
           2232*c[1]**2*c[2]**2 + 252*c[2]**2*c[0]**2 + 1296*c[2]**3*c[0] + \
           3348*c[2]**4 + 24*c[0]**3*c[2] + 3*c[0] - theta[1]
    
    return np.array(f)

def jacobian(c):
    A = [[0]*3 for i in range(3)]
    A[0][0] = 2*c[0] 
    A[0][1] = 4*c[1]
    A[0][2] = 12*c[2]
    A[1][0] = 12*c[1]*c[0] + 36*c[1]*c[2]
    A[1][1] = 24*c[1]**2 + 6*c[0]**2 + 36*c[0]*c[2] + 108*c[2]**2
    A[1][2] = 36*c[0]*c[1] + 216*c[1]*c[2]
    A[2][0] = 120*c[1]**2*c[0] + 576*c[1]**2*c[2] + 504*c[2]**2*c[0] + \
              1296*c[2]**3 + 72*c[0]**2*c[2] + 3
    A[2][1] = 240*c[1]**3 + 120*c[1]*c[0]**2 + 1152*c[1]*c[0]*c[2] + \
              4464*c[1]*c[2]**2
    A[2][2] = 576*c[1]**2*c[0] + 4464*c[1]**2*c[2] + 504*c[2]*c[0]**2 + \
              3888*c[2]**2*c[0] + 13392*c[2]**3 + 24*c[0]**3

    return np.array(A)

def newton(theta, niter = 1e6, tol = 1e-5):
    X=[1,0,0]
    it = 0
    for i in range(int(niter)):
        J = jacobian(X)
        Y = function(X, theta) 
        dX = np.linalg.solve(J,Y) 
        X -= dX 
        it += 1
        if np.linalg.norm(dX) < tol: 
            print('converged')
            break
    print(f"{it=}")
    return X

def broyden(theta, niter = 1e6, tol = 1e-5):
    X0 = [1,0,0]
    B = jacobian(X0)
    it = 0
    for i in range(int(niter)):
        Y = function(X0, theta)
        dX = np.linalg.solve(B,Y)
        X1 = X0 - dX
        Y = function(X1, theta) - function(X0, theta)
        it += 1
        if np.linalg.norm(dX) < tol:
            print('converged')
            break
        else:
            Y = Y.reshape(-1,1)
            dX = dX.reshape(-1,1)
            B = B - (np.matmul((Y + np.matmul(B,dX)), dX.T)) / np.matmul(dX.T,dX)
        X0 = X1
    print(f"{it=}")
    return X1

def main():

    root = Tk()
    root.title("task1")
    w = 30

    i_ICOD = Entry(root, width=w)
    Label(root, text="ICOD : ").grid(row=1, column=0)
    i_ICOD.insert(0, "1")
    i_ICOD.grid(row=1, column=1)

    i_theta1 = Entry(root, width=w)
    Label(root, text="theta1 : ").grid(row=2, column=0)
    i_theta1.insert(0,"0.5")
    i_theta1.grid(row=2, column=1)

    i_theta2 = Entry(root, width=w)
    Label(root, text="theta2 : ").grid(row=3, column=0)
    i_theta2.insert(0,"4.25")
    i_theta2.grid(row=3, column=1)

    i_TOLm = Entry(root, width=w)
    Label(root, text="TOLm : ").grid(row=5, column=0)
    i_TOLm.insert(0, "1e-5")
    i_TOLm.grid(row=5, column=1)

    def myClick():
        global myLabel
        ICOD = json.loads(i_ICOD.get())
        theta = json.loads(i_theta1.get()), json.loads(i_theta2.get())
        TOLm = json.loads(i_TOLm.get())
        if ICOD == 1:
            c = newton(theta, tol = TOLm)
            c = list(map(lambda x: round(x,8), c))
            result = f"c2 = {c[0]}\nc3 = {c[1]}\nc4 = {c[2]}"
        elif ICOD == 2:
            c = broyden(theta, tol = TOLm)
            c = list(map(lambda x: round(x,8), c))
            result = f"c2 = {c[0]}\nc3 = {c[1]}\nc4 = {c[2]}"
        else:
            result = "ICOD não válido"
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

if __name__ == "__main__":
    main()