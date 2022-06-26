import json
import matplotlib.pyplot as plt
from math import sin
from tkinter import Label, Tk, Button, NORMAL, DISABLED, Entry
import pandas as pd

def func(t,x,v,a,w,coef):
    a1,a2,a3 = a
    w1,w2,w3 = w
    m,c,k = coef
    F = a1*sin(w1*t) + a2*sin(w2*t) + a3*sin(w3*t)
    result = (-c*v -k*x + F)/m
    return result

def RKN(h,T,a,w,coef):
    N = T/h
    t_a, x_a, v_a = 0, 0, 0
    t_vec, x_vec, v_vec, f_vec = [], [], [], []
    for k in range(int(N)):
        t = k*h
        k1 = 0.5*h*func(t_a, x_a, v_a, a, w, coef)
        K = 0.5*h*(v_a + 0.5*k1)
        k2 = 0.5*h*func(t_a + 0.5*h, x_a + K, v_a + k1, a, w, coef)
        k3 = 0.5*h*func(t_a + 0.5*h, x_a + K, v_a + k2, a, w, coef)
        L = h*(v_a + k3)
        k4 = 0.5*h*func(t_a + h, x_a + L, v_a + 2*k3, a, w, coef)
        x = x_a + h*(v_a + (1/3)*(k1+k2+k3))
        v = v_a + (1/3)*(k1 + 2*k2 + 2*k3 + k4)
        f = func(t, x, v , a, w, coef)
        f_vec.append(f)
        t_vec.append(round(t,4))
        x_vec.append(x)
        v_vec.append(v)
        t_a, x_a, v_a = t, x, v
    return t_vec, x_vec, v_vec, f_vec

def main():
    root = Tk()
    root.title("task3")
    w = 35
    i_coef = Entry(root, width=w)
    Label(root, text="m, c, k : ").grid(row=1, column=0)
    i_coef.insert(0, "[1.5, 0.1, 2.5]")
    i_coef.grid(row=1, column=1)

    i_a = Entry(root, width=w)
    Label(root, text="a1, a2, a3 : ").grid(row=2, column=0)
    i_a.insert(0,"[-1.87, 1.75, -1.85]")
    i_a.grid(row=2, column=1)

    i_w = Entry(root, width=w)
    Label(root, text="w1, w2, w3 : ").grid(row=3, column=0)
    i_w.insert(0,"[0.075, 1.2, 2.5]")
    i_w.grid(row=3, column=1)

    i_h = Entry(root, width=w)
    Label(root, text="Passo de int : ").grid(row=4, column=0)
    i_h.insert(0, "0.5")
    i_h.grid(row=4, column=1)

    i_T = Entry(root, width=w)
    Label(root, text="Tempo de int : ").grid(row=5, column=0)
    i_T.insert(0, "370")
    i_T.grid(row=5, column=1)

    def myClick():
        global myLabel
        coef = json.loads(i_coef.get())
        a = json.loads(i_a.get())
        w = json.loads(i_w.get())
        h = json.loads(i_h.get())
        T = json.loads(i_T.get())
        result = RKN(h,T,a,w,coef)
        #f = open(f"dados_task3.dat","w+")
        d = {'col1': result[0], 'col2': result[1], 'col3': result[2], 'col4': result[3]}
        df = pd.DataFrame(data=d)
        df.to_csv(f"dados_task3.txt", header=None, index=None, sep='\t', mode='a')
        fig, axs = plt.subplots(3, sharex=True, sharey=True)
        axs[0].plot(result[0],result[1])
        axs[0].set_title("Deslocamento")
        axs[1].plot(result[0],result[2])
        axs[1].set_title("Velocidade")
        axs[2].plot(result[0],result[3])
        axs[2].set_title("Aceleração")
        axs[0].set(xlabel='t', ylabel='y(t)')
        axs[1].set(xlabel='t', ylabel='y\'(t)')
        axs[2].set(xlabel='t', ylabel='y\'\'(t)')
        plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.4, 
                        hspace=0.4)
        plt.ion()
        plt.show()
        myButton['state'] = DISABLED

    def Deletar():
        plt.close()
        myButton['state'] = NORMAL

    myButton = Button(root, text = "Calcular", command = myClick)
    myButton.grid(row=6, column=1)

    myButton2 = Button(root, text = "Resetar", command = Deletar)
    myButton2.grid(row=7, column=1)

    root.mainloop()

if __name__ == "__main__":
    main()