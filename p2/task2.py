import numpy as np
import math
import json
from tkinter import Label, Tk, Button, NORMAL, DISABLED, Entry

W_Legendre =      {2:  {"p": [-0.5773502691896257, 0.5773502691896257], "w": [1.0, 1.0]},
                   3:  {"p": [0.0, -0.7745966692414834, 0.7745966692414834], "w": [0.8888888888888888, 0.5555555555555556, 0.5555555555555556]},
                   4:  {"p": [-0.3399810435848563, 0.3399810435848563, -0.8611363115940526, 0.8611363115940526], "w": [0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538]},
                   5:  {"p": [0.0, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640, 0.9061798459386640], "w": [0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891]},
                   6:  {"p": [0.6612093864662645, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, -0.9324695142031521, 0.9324695142031521], "w": [0.3607615730481386, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.1713244923791704, 0.1713244923791704]},
                   7:  {"p": [0.0, 0.4058451513773972, -0.4058451513773972, -0.7415311855993945, 0.7415311855993945, -0.9491079123427585, 0.9491079123427585], "w": [0.4179591836734694, 0.3818300505051189, 0.3818300505051189, 0.2797053914892766, 0.2797053914892766, 0.1294849661688697, 0.1294849661688697]},
                   8:  {"p": [-0.1834346424956498, 0.1834346424956498, -0.5255324099163290, 0.5255324099163290, -0.7966664774136267, 0.7966664774136267, -0.9602898564975363, 0.9602898564975363], "w": [0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.3137066458778873, 0.2223810344533745, 0.2223810344533745, 0.1012285362903763, 0.1012285362903763]},
                   9:  {"p": [0.0, -0.8360311073266358, 0.8360311073266358, -0.9681602395076261, 0.9681602395076261, -0.3242534234038089, 0.3242534234038089, -0.6133714327005904, 0.6133714327005904], "w": [0.3302393550012598, 0.1806481606948574, 0.1806481606948574, 0.0812743883615744, 0.0812743883615744, 0.3123470770400029, 0.3123470770400029, 0.2606106964029354, 0.2606106964029354]},
                   10: {"p": [-0.1488743389816312, 0.1488743389816312, -0.4333953941292472, 0.4333953941292472, -0.6794095682990244, 0.6794095682990244, -0.8650633666889845, 0.8650633666889845, -0.9739065285171717, 0.9739065285171717], "w": [0.2955242247147529, 0.2955242247147529, 0.2692667193099963, 0.2692667193099963, 0.2190863625159820, 0.2190863625159820, 0.1494513491505806, 0.1494513491505806, 0.0666713443086881, 0.0666713443086881]}}

def func(x,coef):
    c1, c2, c3, c4 = coef
    return c1*math.e**(c2*x) + c3*x**c4

def bissecao(coef, intervalo, tol = 1e-5):
    a, b = intervalo
    aux1 = 0
    aux2 = 0
    it = 0
    while abs(b-a) > tol:
        x = (a+b) / 2.0
        f = func(x,coef)
        aux1 +=1
        it +=1
        if f > 0.0:
            b = x
            aux2 +=1
        else:
            a = x
            aux2 -=1
    if aux1 == abs(aux2):
        return "[Raiz(es) fora do intervalo]"
    print(f"{it=}")
    return round(x, 8)

def newton(coef, intervalo, niter=1e6, tol=1e-5):
    a, b = intervalo
    def func_prime(a,coef):
        h = 1e-8
        return (func(a+h,coef) - func(a,coef)) / h
    x0 = (a+b)/2
    it = 0
    for _ in range(int(niter)):
        x1 = x0 - func(x0, coef) / func_prime(x0, coef)
        tolk = abs(x1-x0)
        it += 1
        if tolk < tol:
            print(f"{it=}")
            return round(x1, 8)
        x0 = x1
    return "[Não convergiu]"

def quadratura_polinomial(coef, intervalo, num_int):
    a, b = intervalo
    delta_x = float(abs(b-a)) / (num_int-1)
    b_vec = np.array([(b**i-a**i)/i for i in range(1, num_int+1)])
    x_vec = np.array([a+(i-1)*delta_x for i in range(1, num_int+1)])
    vandermond_matrix = np.array([[x_vec[j]**i for j in range(num_int)] for i in range(num_int)])
    w = np.linalg.solve(vandermond_matrix, b_vec)
    result = 0
    for i in range(num_int):
        result += w[i]*func(x_vec[i], coef)
    print(w)
    return round(result, 8)

def gauss_legendre(coef, intervalo, num_int):
    a, b = intervalo
    w_legendre = W_Legendre.get(num_int).get("w")
    p_legendre = W_Legendre.get(num_int).get("p")
    L = b - a
    x_vec = np.array([(a+b+p_legendre[i]*L)/2 for i in range(num_int)])
    integral_sum = 0
    for i in range(num_int):
        integral_sum += func(x_vec[i], coef)*w_legendre[i]
    result = integral_sum*L / 2
    return round(result, 8)

def diferencas_finitas(coef, a, delta_x, num):
    if num == 1:
        result = (func(a+delta_x, coef)-func(a, coef)) / delta_x
        return round(result, 8)
    elif num == 2:
        result = (func(a, coef)-func(a-delta_x, coef)) / delta_x
        return round(result, 8)
    elif num == 3:
        result = (func(a+delta_x, coef) - func(a-delta_x, coef)) / (2*delta_x)
        return round(result, 8)
    return "Método inválido"

def richard(coef, a, delta_x1, delta_x2):
    f1 = diferencas_finitas(coef, a, delta_x1, 1)
    f2 = diferencas_finitas(coef, a, delta_x2, 1)
    return round(2*f2 - f1, 8)

def main():
    root = Tk()
    root.title("task2")
    root.geometry("310x50+400+300")

    def f_icod1():
        global i11,i12,i13,i14,ICOD
        ICOD = 1
        destroy_init()
        Label(root, text="Tarefa escolhida: 1 - Raiz").grid(row=0, column=1)
        i11 = insert_input("Constantes: ", 1, 0, "[0.25, 1, -2, 1.3]")
        i12 = insert_input("ICOD: ", 2, 0, "2")
        i13 = insert_input("Intervalo: ", 3, 0, "[1, 5]")
        i14 = insert_input("TOLm: ", 4, 0, "1e-5")
        b_calcular.grid(row=5, column=1)
        b_resetar.grid(row=7, column=1)
        b_alterar.grid(row=8,column=1)
        root.geometry("310x230")

    def f_icod2():
        global i21,i22,i23,i24,ICOD
        ICOD = 2
        destroy_init()
        Label(root, text="Tarefa escolhida: 2 - Integral").grid(row=0, column=1)
        i21 = insert_input("Constantes: ", 1, 0, "[0.25, 1, -2, 1.3]")
        i22 = insert_input("ICOD: ", 2, 0, "2")
        i23 = insert_input("Intervalo: ", 3, 0, "[0, 4]")
        i24 = insert_input("Pontos de Int: ", 4, 0, "10")
        b_calcular.grid(row=5, column=1)
        b_resetar.grid(row=7, column=1)
        b_alterar.grid(row=8,column=1)
        root.geometry("310x230")

    def f_icod3():
        global i31,i32,i33,i34,ICOD
        ICOD = 3
        destroy_init()
        Label(root, text="Tarefa escolhida: 3 - Derivada").grid(row=0, column=1)
        i31 = insert_input("Constantes: ", 1, 0, "[0.25, 1, -2, 1.3]")
        i32 = insert_input("ICOD: ", 2, 0, "2")
        i33 = insert_input("Ponto: ", 3, 0, "10")
        i34 = insert_input("delta_x: ", 4, 0, "1e-8")
        b_calcular.grid(row=5, column=1)
        b_resetar.grid(row=7, column=1)
        b_alterar.grid(row=8,column=1)
        root.geometry("310x230")

    def f_icod4():
        global i41,i42,i43,i44,ICOD
        ICOD = 4
        destroy_init()
        Label(root, text="Tarefa escolhida: 4 - Derivada RE").grid(row=0, column=1)
        i41 = insert_input("Constantes: ", 1, 0, "[0.25, 1, -2, 1.3]")
        i42 = insert_input("Ponto: ", 2, 0, "10")
        i43 = insert_input("delta_x1: ", 3, 0, "1e-6")
        i44 = insert_input("delta_x2: ", 4, 0, "1e-8")
        b_calcular.grid(row=5, column=1)
        b_resetar.grid(row=7, column=1)
        b_alterar.grid(row=8,column=1)
        root.geometry("310x230")

    b_icod1 = Button(root, text = "1", width=3,command=f_icod1)
    b_icod2 = Button(root, text = "2", width=3,command=f_icod2)
    b_icod4 = Button(root, text = "4", width=3,command=f_icod4)
    b_icod3 = Button(root, text = "3", width=3,command=f_icod3)
    init_label = Label(root, text="Tarefa (ICOD): ")

    def f_calcular():
        global l_result
        if ICOD == 1:
            coef = json.loads(i11.get())
            sub_ICOD = json.loads(i12.get())
            intervalo = json.loads(i13.get())
            tol = json.loads(i14.get())
            if sub_ICOD == 1:
                result = f"Raiz = {bissecao(coef, intervalo, tol=tol)}"
            elif sub_ICOD == 2:
                result = f"Raiz = {newton(coef, intervalo, tol=tol)}"
            else:
                result = "Método inválido"
        elif ICOD == 2:
            coef = json.loads(i21.get())
            sub_ICOD = json.loads(i22.get())
            intervalo = json.loads(i23.get())
            num_int = json.loads(i24.get())
            if sub_ICOD == 1:
                result = f"Integral em {intervalo} = {quadratura_polinomial(coef, intervalo, num_int)}"
            elif sub_ICOD == 2:
                result = f"Integral em {intervalo} = {gauss_legendre(coef, intervalo, num_int)}"
            else:
                result = "Método inválido"
        elif ICOD == 3:
            coef = json.loads(i31.get())
            sub_ICOD = json.loads(i32.get())
            a = json.loads(i33.get())
            delta_x = json.loads(i34.get())
            if 1 <= sub_ICOD <= 3:
                result = f"Derivada em {a} = {diferencas_finitas(coef,a,delta_x,sub_ICOD)}"
            else:
                result = "Método inválido"
        elif ICOD == 4:
            coef = json.loads(i41.get())
            a = json.loads(i42.get())
            delta_x1 = json.loads(i43.get())
            delta_x2 = json.loads(i44.get())
            result = f"Derivada em {a} = {richard(coef, a, delta_x1, delta_x2)}" 
        l_result = Label(root, text = result)
        l_result.grid(row = 6, column = 1)
        b_calcular['state'] = DISABLED

    def f_resetar():
        l_result.destroy()
        b_calcular['state'] = NORMAL

    b_calcular = Button(root, text = "Calcular", command=f_calcular)
    b_resetar = Button(root, text = "Resetar",command=f_resetar)
    b_alterar = Button(root, text = "Alterar tarefa", command=lambda:[root.destroy(),main()])
    
    init_label.grid(row=1, column=0)
    b_icod1.grid(row=1, column=1)
    b_icod2.grid(row=1, column=2)
    b_icod3.grid(row=1, column=3)
    b_icod4.grid(row=1, column=4)

    def destroy_init():
        b_icod1.destroy()
        b_icod2.destroy()
        b_icod3.destroy()
        b_icod4.destroy()
        init_label.destroy()

    def insert_input(i_text, r, c, i_insert):
        entry = Entry(root, width=37)
        label = Label(root, text=i_text)
        entry.insert(0, i_insert)
        label.grid(row=r, column=c)
        entry.grid(row=r, column=str(int(c+1)))
        return entry

    root.mainloop()

if __name__ == "__main__":
    main()