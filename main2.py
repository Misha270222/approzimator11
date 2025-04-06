import math
import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from scipy.optimize import root, root_scalar
import traceback


class Tab2:
    def __init__(self, master, data_mu1, data_mu2, main_entries):
        self.tab = ttk.Frame(master)
        if not isinstance(data_mu1, dict) or not isinstance(data_mu2, dict):
            raise ValueError("Некорректные данные для μ1 и μ2!")
        self.data_mu1 = data_mu1
        self.data_mu2 = data_mu2
        self.main_entries = main_entries
        self.a = 0
        self.b = 0
        self.gamma = (self.b - self.a)/2
        self.create_widgets()

    def create_widgets(self):
        # Поля ввода
        fields = [
            ("Внутренний радиус a", "a"),
            ("Внешний радиус b", "b"),
            ("Нагрузка q_a", "qa"),
            ("Нагрузка q_b", "qb"),
            ("Коэф. Пуассона ν", "nu"),
            ("Эпсилон", "epsilon0"),
        ]

        self.entries = {}
        for i, (field, name) in enumerate(fields):
            tk.Label(self.tab, text=field).grid(row=i, column=0)
            entry = tk.Entry(self.tab)
            entry.grid(row=i, column=1)
            self.entries[name] = entry

        # Установка начальных значений
        self.entries["a"].insert(0, "1")
        self.entries["b"].insert(0, "3")
        self.entries["qa"].insert(0, "0.001")
        self.entries["qb"].insert(0, "0.03")
        self.entries["nu"].insert(0, "0.4")
        self.entries["epsilon0"].insert(0, "0.03")

        # Дополнительные поля для границы слоев
        # tk.Label(self.tab, text="Гамма").grid(row=len(fields), column=0)
        self.entry_boundary = tk.Entry(self.tab)
        self.entry_boundary.grid(row=len(fields), column=1)
        self.entry_boundary.insert(0, "2")

        # Графики
        self.fig = plt.figure(figsize=(14, 8))
        gs = self.fig.add_gridspec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])

        # Первый столбец для первых двух графиков
        self.ax_u = self.fig.add_subplot(gs[0, 0])  # Верхний график
        self.ax_sr = self.fig.add_subplot(gs[1, 0])  # Нижний график

        # Второй столбец для третьего графика (занимает обе строки)
        self.ax_st = self.fig.add_subplot(gs[0, 1])
        self.ax_stheta = self.fig.add_subplot(gs[1, 1])  # Полный столбец

        # Отображение графиков в Tkinter
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.tab)
        self.canvas.get_tk_widget().grid(row=len(fields) + 2, columnspan=2)

        # Настройка макета
        self.fig.tight_layout()

        # Кнопка расчета
        btn_calculate = tk.Button(self.tab, text="Рассчитать", command=self.update_results)
        btn_calculate.grid(row=len(fields) + 1, columnspan=2)

    def get_current_mu_func(self):
        current_mu = self.main_entries['var_mu'].get()
        if current_mu == 1:
            return self.data_mu1["function"]
        return self.data_mu2["function"]

    def mu1(self, r):
        if not self.data_mu1["function"]:
            raise ValueError("Функция μ1 не определена!")
        # print("mu1: ", self.data_mu1)
        return self.data_mu1["function"](r)

    def mu2(self, r):
        if not self.data_mu2["function"]:  # Исправлено на mu2
            raise ValueError("Функция μ2 не определена!")
        return self.data_mu2["function"](r)

    def X(self, gamma):
        a_val = float(self.entries["a"].get())  # Получаем актуальное значение a
        steps = 10000
        t = np.linspace(a_val , gamma, steps + 1)
        x_integrand = t / self.mu1(t)
        h = (gamma - a_val) / steps
        integral = h * (0.5 * (x_integrand[0] + x_integrand[-1]) + np.sum(x_integrand[1:-1]))
        return integral

    def Y(self, gamma):
        b_val = float(self.entries["b"].get())
        steps = 10000
        t = np.linspace(gamma, b_val, steps + 1)
        integrand = (self.mu1(t) + self.mu2(t)) / t ** 3
        h = (b_val - gamma) / steps
        integral = h * (0.5 * (integrand[0] + integrand[-1]) + np.sum(integrand[1:-1]))
        return integral

    def Z(self, gamma):
        b_val = float(self.entries["b"].get())
        steps = 10000
        t = np.linspace(gamma, b_val, steps + 1)
        integrand = self.mu2(t) / t
        h = (b_val - gamma) / steps
        integral = h * (0.5 * (integrand[0] + integrand[-1]) + np.sum(integrand[1:-1]))
        return integral


    def W(self, gamma):
        """Функция W, определенная формулой"""
        y_val = self.Y(gamma)
        z_val = self.Z(gamma)
        mu1_val = self.mu1(gamma)
        return 2 * y_val - ((mu1_val + 2 * z_val) / (gamma**2))

    def Omega1(self, gamma):
        """Компонента Ω1 для итерационного процесса"""
        # print("o1 begin")
        y_val = self.Y(gamma)
        # print("y o1 ok")
        w_val = self.W(gamma)

        qa = float(self.entries["qa"].get())
        qb = float(self.entries["qb"].get())
        a_val = float(self.entries["a"].get())
        mu1_val = self.mu1(a_val)
        b_val = float(self.entries["b"].get())
        epsilon0 = float(self.entries["epsilon0"].get())

        return y_val*epsilon0*((gamma)**2) + ((qa*(a_val**2)) / (2*mu1_val)) * w_val + qb/2

    def Omega2(self, gamma):
        """Компонента Ω2 для итерационного процесса"""
        x_val = self.X(gamma)
        w_val = self.W(gamma)
        z_val = self.Z(gamma)

        nu = float(self.entries["nu"].get())
        a_val = float(self.entries["a"].get())
        mu1_a_val = self.mu1(a_val)
        mu1_gamma_val = self.mu1(gamma)

        # p1 = (1 - nu) / (1 - 2 * nu)
        # p2 = (x_val + ((1 - nu) / (1 - 2 * nu)) * (a_val**2 / mu_val))
        return (1 - nu) / (1 - 2 * nu) + (x_val + ((1 - nu) / (1 - 2 * nu)) * ((a_val**2) / mu1_a_val))*w_val + z_val/mu1_gamma_val

    def Omega(self, gamma):
        """Объединенная функция Ω для определения γ"""
        # print("omega begin")
        omega1 = self.Omega1(gamma)
        # print("omega1 ok")
        mu_val = self.mu1(gamma)
        # print("mu ok")
        omega2 = self.Omega2(gamma)
        # print("omega2 ok")
        epsilon0 = float(self.entries["epsilon0"].get())
        return omega1 - epsilon0*mu_val*omega2

    # def gamma_equation(self, gamma):
    #     """Уравнение для нахождения γ методом root"""
    #     return self.Omega(self.r_boundary, gamma)

    def ef(self, gamma):
        X_g = self.X(gamma)
        return 2/(gamma**2) * X_g - 1/(self.mu1(gamma))
    def find_gamma(self):
        try:
            a = float(self.entries['a'].get())
            b_val = float(self.entries['b'].get())
            # Используем метод Брента для быстрого нахождения корня
            result = root_scalar(self.Omega, bracket=[a, b_val], method='brentq', xtol=1e-6)
            if not result.converged:
                raise RuntimeError("Решение не найдено.")
            self.gamma = result.root
            return self.gamma
        except ValueError as e:
            raise RuntimeError(f"Ошибка: {str(e)}") from e

    def find_integration_constants(self):
        a = float(self.entries["a"].get())
        b = float(self.entries["b"].get())
        qa = float(self.entries["qa"].get())
        qb = float(self.entries["qb"].get())
        v = float(self.entries["nu"].get())
        epsilon0 = float(self.entries["epsilon0"].get())

        X_gamma = self.X(self.gamma)
        Y_gamma = self.Y(self.gamma)
        Z_gamma = self.Z(self.gamma)

        mu1_gamma = self.mu1(self.gamma)

        phi_gamma11 = X_gamma
        phi_gamma21 = ((1 - v)/(1 - 2*v))*((a**2) / (self.mu1(a)))
        # phi_gamma31 = ((1 - v) / (1 - 2*v)) - (X_gamma / (self.gamma)**2) * (mu1_gamma + 2 * Z_gamma) + ((Z_gamma)/(mu1_gamma))
        # phi_gamma32 = (-1/((self.gamma)**2))*(mu1_gamma + 2 * Z_gamma)
        # phi_gamma33 = 2 * Y_gamma

        d_gamma1 = (-1/2)*epsilon0*((self.gamma)**2)
        d_gamma2 = ((-1*qa)/(2*self.mu1(a)))*(a**2)

        # d2 = (-1*qb)/2

        W_gamma = self.W(self.gamma)
        Omega1_gamma = self.Omega1(self.gamma)
        Omega2_gamma = self.Omega2(self.gamma)
        Omega_gamma = self.Omega(self.gamma)

        A1 = -1 * ((Omega1_gamma)/(Omega2_gamma))
        B1 = A1 * phi_gamma21 - d_gamma2
        C1 = A1 * phi_gamma11 + B1 - d_gamma1

        # print(f"A1 = {A1}, B1 = {B1}, C1 = {C1}")
        # print(f"a = {a}, b = {b}, qa = {qa}, qb = {qb}, nu = {v}, epsiloin0 = {epsilon0}, gamma = {self.gamma}")
        # print(f"X_g = {X_gamma}")
        # print(f"Y_g = {Y_gamma}")
        # print(f"Z_g = {Z_gamma}")
        # print(f"mu1_g = {mu1_gamma}")
        # print(f"Omega1_g = {Omega1_gamma}")
        # print(f"Omega2_g = {Omega2_gamma}")
        return A1, B1, C1

    def calculate_displacements(self):
        self.find_gamma()
        A1, B1, C1 = self.find_integration_constants()
        epsilon0 = float(self.entries["epsilon0"].get())
        a_rad = float(self.entries["a"].get())
        b_rad = float(self.entries["b"].get())
        r = np.linspace(a_rad, b_rad, 100)
        v = float(self.entries["nu"].get())
        qb = float(self.entries["qb"].get())

        # Векторизованные вычисления
        mask_inner = r <= self.gamma
        r_inner = r[mask_inner]
        r_outer = r[~mask_inner]

        # Перемещения
        u = np.zeros_like(r)
        u[mask_inner] = (A1 * np.array([self.X(ri) for ri in r_inner]) + B1) / r_inner
        u[~mask_inner] = C1 / r_outer - 0.5 * epsilon0 * r_outer

        # Деформации
        theta = np.zeros_like(r)
        theta[mask_inner] = A1 / self.mu1(r_inner)

        # Напряжения
        sigma = np.zeros_like(r)
        sigma[mask_inner] = 2 * self.mu1(r_inner) * (
            A1 * ((1 - v) / (1 - 2 * v) / self.mu1(r_inner) - np.array([self.X(ri) for ri in r_inner]) / r_inner**2) -
            B1 / r_inner**2
        )

        X_g = self.X(self.gamma)
        term = 2 * (A1 * (2 * X_g / self.gamma ** 2 - 1 / self.mu1(self.gamma)) + 2 * B1 / self.gamma ** 2)
        sigma[~mask_inner] = term * self.Z(r_outer) - 4 * C1 * self.Y(r_outer) - float(self.entries["qb"].get())

        sigma_theta = np.zeros_like(r)
        if len(r_inner) > 0:
            mu1_inner = self.mu1(r_inner)
            X_inner = np.array([self.X(ri) for ri in r_inner])
            sigma_theta[mask_inner] = 2 * mu1_inner * (
                A1 * (v/1-2*v)/mu1_inner +
                (X_inner + B1)/(r_inner**2)
            )

        if len(r_outer) > 0:
            mu_outer = self.mu1(r_outer) + self.mu2(r_outer)
            Y_outer = np.array([self.Y(ro) for ro in r_outer])
            Z_outer = np.array([self.Z(ro) for ro in r_outer])
            mu0_outer = self.mu2(r_outer)
            ef_g = self.ef(self.gamma)

            sigma_theta[~mask_inner] = (
                4 * C1 * (mu_outer/r_outer**2 - Y_outer) +
                2 * (A1 * ef_g + 2*B1/self.gamma**2) + (Z_outer - mu0_outer) -
                qb
            )
        # Для внешней части вычисляем один раз общие параметры

        return r, u, theta, sigma, sigma_theta
        # return r, u, sigma_r, sigma_theta

    def update_results(self):
        try:
            # ИЗМЕНЕНИЕ 7: Добавить проверку функций
            if not self.data_mu1["function"] or not self.data_mu2["function"]:
                raise ValueError("Сначала выполните аппроксимацию μ1 и μ2 в первой вкладке!")

            self.r_boundary = float(self.entry_boundary.get())
            r, u, theta, sigma, sigma_theta = self.calculate_displacements()

            # Граница между слоями для отображения на графиках
            boundary = float(self.gamma)

            # Очистка графиков
            self.ax_u.clear()
            self.ax_sr.clear()
            self.ax_st.clear()
            self.ax_stheta.clear()

            # Построение графиков с выделением границы слоев
            self.ax_u.plot(r, u)
            self.ax_u.axvline(x=boundary, color='r', linestyle='--', label='Граница слоев')
            self.ax_u.set_title("перемещение u(r)")
            self.ax_u.set_xlabel("Радиус r")
            self.ax_u.set_ylabel("Перемещение")
            self.ax_u.grid(True)
            self.ax_u.legend()

            self.ax_sr.plot(r, theta)
            eps = float(self.entries["epsilon0"].get())
            eps *= -1
            self.ax_sr.axhline(y=eps, color='r', linestyle='--', label='ε₀')
            self.ax_sr.set_title("Деформации Θ(r)")
            self.ax_sr.set_xlabel("Радиус r")
            self.ax_sr.set_ylabel("Деформации")
            self.ax_sr.grid(True)
            self.ax_sr.legend()

            self.ax_st.plot(r, sigma)
            self.ax_st.axvline(x=boundary, color='r', linestyle='--', label='Граница слоев')
            self.ax_st.set_title("напряжение σ(r)")
            self.ax_st.set_xlabel("Радиус r")
            self.ax_st.set_ylabel("Напряжение")
            self.ax_st.grid(True)
            self.ax_st.legend()

            self.ax_stheta.plot(r, sigma_theta)
            self.ax_stheta.axvline(x=boundary, color='r', linestyle='--', label='Граница слоев')
            self.ax_stheta.set_title("напряжение σ_θ(r)")
            self.ax_stheta.set_xlabel("Радиус r")
            self.ax_stheta.set_ylabel("Напряжение")
            self.ax_stheta.grid(True)
            self.ax_stheta.legend()

            # Обновление графиков
            self.fig.tight_layout()
            self.canvas.draw()

            # Вывод диагностической информации
            print(f"Радиус: от {r[0]} до {r[-1]}")
            print(f"Перемещение: от {min(u)} до {max(u)}")
            # print(f"Радиальное напряжение: от {min(sigma_r)} до {max(sigma_r)}")
            # print(f"Окружное напряжение: от {min(sigma_theta)} до {max(sigma_theta)}")


        except Exception as e:

            # ИЗМЕНЕНИЕ 8: Показывать ошибку в GUI

            import tkinter.messagebox as mb

            mb.showerror("Ошибка расчета", f"{str(e)}")

            print(f"Ошибка: {traceback.format_exc()}")