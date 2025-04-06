import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from scipy.optimize import curve_fit
from main2 import Tab2
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk
)
from style import configure_styles

global zoom_active, pan_active, prev_x, prev_y
zoom_active = False
pan_active = False
prev_x = None
prev_y = None

entry_a = None
entry_b = None
# Глобальные переменные для хранения данных μ1 и μ2
data_mu1 = {
    "mu0": "",
    "function": None,  # Здесь будет lambda функция
    "params": [],      # Параметры аппроксимации
    "table_data": [],
}
data_mu2 = {
    "mu0": "",
    "function": None,  # Здесь будет lambda функция
    "params": [],      # Параметры аппроксимации
    "table_data": [],
}

# Список для хранения всех линий графиков
all_lines = []

# Переменная для управления отображением точек
show_points = True

def on_mousewheel(event):
    # Зуммирование по прокрутке мыши
    base_scale = 1.2  # Коэффициент масштабирования
    if event.button == 'up':
        scale_factor = 1 / base_scale
    else:
        scale_factor = base_scale

    # Получаем текущие пределы осей
    current_xlim = ax.get_xlim()
    current_ylim = ax.get_ylim()

    # Вычисляем новую ширину и высоту осей
    new_width = (current_xlim[1] - current_xlim[0]) * scale_factor
    new_height = (current_ylim[1] - current_ylim[0]) * scale_factor

    # Центрируем масштабирование относительно позиции курсора
    x_center = event.xdata
    y_center = event.ydata
    if x_center is not None and y_center is not None:
        ax.set_xlim([x_center - new_width/2, x_center + new_width/2])
        ax.set_ylim([y_center - new_height/2, y_center + new_height/2])
        canvas.draw()


def on_press(event):
    # Начало панорамирования при нажатии левой кнопки мыши
    global pan_active, prev_x, prev_y
    if event.button == 1:  # Левая кнопка
        pan_active = True
        prev_x = event.xdata
        prev_y = event.ydata

def on_release(event):
    # Окончание панорамирования
    global pan_active
    pan_active = False

SMOOTHING_FACTOR = 0.5  # Значение от 0 до 1 (0.5 — усреднение на 50%)
def on_motion(event):
    global pan_active, prev_x, prev_y
    if pan_active and prev_x is not None and prev_y is not None:
        if event.xdata is None or event.ydata is None:
            return
        dx = (event.xdata - prev_x) * SMOOTHING_FACTOR
        dy = (event.ydata - prev_y) * SMOOTHING_FACTOR
        current_xlim = ax.get_xlim()
        current_ylim = ax.get_ylim()
        ax.set_xlim(current_xlim[0] - dx, current_xlim[1] - dx)
        ax.set_ylim(current_ylim[0] - dy, current_ylim[1] - dy)
        prev_x = event.xdata
        prev_y = event.ydata
        canvas.draw()

def save_data():
    data = {
        "μ0": entries['mu0'].get(),
        "a": entry_a.get(),
        "b": entry_b.get(),
        "ν1": entry_nu1.get(),
        "ν2": entry_nu2.get(),
        "qa": entry_qa.get(),
        "qb": entry_qb.get(),
        "ε0": entry_epsilon0.get(),
        "μ1": combo_mu1.get() if var_mu.get() == 1 else None,
        "μ2": combo_mu2.get() if var_mu.get() == 2 else None,
    }
    with open("output.txt", "w", encoding="utf-8") as file:
        for key, value in data.items():
            file.write(f"{key}: {value}\n")
    messagebox.showinfo("Сохранение", "Данные успешно сохранены в файл output.txt")

def update_plot(auto=False):
    global show_points
    r_values = []
    mu_values = []
    # Считываем данные из таблицы
    for row in table.get_children():
        values = table.item(row, "values")
        try:
            r = float(values[0])
            mu = float(values[1])
            r_values.append(r)
            mu_values.append(mu)
        except (ValueError, IndexError):
            messagebox.showerror("Ошибка", "Некорректные данные в таблице.")
            return
    r_values = np.array(r_values)
    mu_values = np.array(mu_values)
    if len(r_values) < 2 or len(mu_values) < 2:
        messagebox.showerror("Ошибка", "Необходимо ввести хотя бы две точки данных.")
        return

    selected_func = combo_mu1.get() if var_mu.get() == 1 else combo_mu2.get()

    # Определяем функцию для аппроксимации
    def get_function(selected_func):
        a = float(entry_a.get()) if entry_a.get() else 0
        b = float(entry_b.get()) if entry_b.get() else 0

        if selected_func == "μ(r) = μ₀(1 - k _exp(-α_ r))":
            def func(r, k, alpha):
                mu0 = float(entries['mu0'].get())
                return mu0 * (1 - k * np.exp(-alpha * r))
            return func, [1, 1]

        elif selected_func == "μ(r) = μ₀(1 + k_ ((b - r) / (b - a))^n)":
            def func(r, k, n):
                mu0 = float(entries['mu0'].get())
                a_val = a
                b_val = b
                return mu0 * (1 + k * ((b_val - r) / (b_val - a_val)) ** n)
            return func, [0.5, 4] # Начальные параметры

        elif selected_func == "μ(r) = μ₀[1 + (k-1)(a/r)^n]":
            def func(r, k, n):
                mu0 = float(entries['mu0'].get())
                a_val = a
                return mu0 * (1 + (k - 1) * (a_val / r) ** n)
            return func, [0.5, 4]  # Начальные параметры

        elif selected_func == "μ(r) = m(1 + k e^{-ar} (cos(βr + ω) + sin(βr + ω)))":
            def func(r, k, beta, omega):
                mu0 = float(entries['mu0'].get())
                a_val = a
                return mu0 * (1 + k * np.exp(-a_val * r) * (np.cos(beta * r + omega) + np.sin(beta * r + omega)))
            return func, [0.1, 2.0, 0.0]  # Начальные параметры

    func, initial_params = get_function(selected_func)

    try:
        popt, pcov = curve_fit(func, r_values, mu_values, p0=initial_params)
        mu_fit = func(r_values, *popt)

        param_names_dict = {
            "μ(r) = μ₀(1 - k exp(-α r))": ["k", "alpha"],
            "μ(r) = μ₀(1 + k_ ((b - r) / (b - a))^n)": ["k", "n", "c", "a"],
            "μ(r) = μ₀[1 + (k-1)(a/r)^n]": ["k", "n", "a"],
            "μ(r) = m(1 + k e^{-ar} (cos(βr + ω) + sin(βr + ω)))": ["k", "a", "beta", "omega"]
        }
        param_names = param_names_dict.get(selected_func, [])
        if param_names:
            print(f"\nАппроксимация функцией: {selected_func}")
            print("Оптимизированные коэффициенты:")
            for name, value in zip(param_names, popt):
                print(f"  {name} = {value:.4f}")
            r_test = 1
            mu_at_dot = func(r_test, *popt)
            print(f"Значение функции в точке {r_test} = {mu_at_dot:.4f}")
        else:
            print(f"Неизвестная функция: {selected_func}")

        # Добавляем новый график
        line, = ax.plot(r_values, mu_fit, label=f'Аппроксимация ({selected_func})')
        all_lines.append(line)

        # Добавляем точки, если show_points == True
        if show_points:
            ax.scatter(r_values, mu_values, color='blue', label='Исходные данные')

        # Настройка осей
        if auto:
            xmin, xmax = min(r_values), max(r_values)
            ymin, ymax = min(mu_values), max(mu_values)
            xstep = (xmax - xmin) / 10
            ystep = (ymax - ymin) / 10
            entry_xmin.delete(0, tk.END)
            entry_xmin.insert(0, f"{xmin:.2f}")
            entry_xmax.delete(0, tk.END)
            entry_xmax.insert(0, f"{xmax:.2f}")
            entry_ymin.delete(0, tk.END)
            entry_ymin.insert(0, f"{ymin:.2f}")
            entry_ymax.delete(0, tk.END)
            entry_ymax.insert(0, f"{ymax:.2f}")
            entry_xstep.delete(0, tk.END)
            entry_xstep.insert(0, f"{xstep:.2f}")
            entry_ystep.delete(0, tk.END)
            entry_ystep.insert(0, f"{ystep:.2f}")
        else:
            xmin = float(entry_xmin.get())
            xmax = float(entry_xmax.get())
            ymin = float(entry_ymin.get())
            ymax = float(entry_ymax.get())
            xstep = float(entry_xstep.get())
            ystep = float(entry_ystep.get())

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.set_xticks(np.arange(xmin, xmax + xstep, xstep))
        ax.set_yticks(np.arange(ymin, ymax + ystep, ystep))
        ax.set_xlabel("r")
        ax.set_ylabel("μ")
        ax.set_title("График зависимости μ от r")
        ax.legend()
        canvas.draw()

    except RuntimeError as e:
        messagebox.showerror("Ошибка", f"Не удалось построить график: {str(e)}")
    else:
        current_mu = var_mu.get()
        mu_data = data_mu1 if current_mu == 1 else data_mu2
        mu_data["function"] = lambda r, f=func, p=popt: f(r, *p)
        mu_data["params"] = popt.tolist()
def switch_mu():
    current_mu = var_mu.get()
    if current_mu == 1:
        combo_mu1.grid(row=len(fields_tab1) + 2, column=0, columnspan=2, pady=5)
        combo_mu2.grid_forget()
        entries['mu0'].delete(0, tk.END)
        entries['mu0'].insert(0, data_mu1["mu0"])
        combo_mu1.set(data_mu1["function"])
        load_table_data(data_mu1["table_data"])
    elif current_mu == 2:
        combo_mu2.grid(row=len(fields_tab1) + 3, column=0, columnspan=2, pady=5)
        combo_mu1.grid_forget()
        entries['mu0'].delete(0, tk.END)
        entries['mu0'].insert(0, data_mu2["mu0"])
        combo_mu2.set(data_mu2["function"])
        load_table_data(data_mu2["table_data"])

def load_table_data(data):
    for i, row in enumerate(table.get_children()):
        if i < len(data):
            table.item(row, values=data[i])
        else:
            table.item(row, values=(0, 0))

def save_current_mu_data():
    current_mu = var_mu.get()
    if current_mu == 1:
        data_mu1["mu0"] = entries['mu0'].get()
        data_mu1["function"] = combo_mu1.get()
        data_mu1["table_data"] = [table.item(row, "values") for row in table.get_children()]
    elif current_mu == 2:
        data_mu2["mu0"] = entries['mu0'].get()
        data_mu2["function"] = combo_mu2.get()
        data_mu2["table_data"] = [table.item(row, "values") for row in table.get_children()]

def clear_plot():
    ax.clear()
    all_lines.clear()
    ax.set_xlabel("r")
    ax.set_ylabel("μ")
    ax.set_title("График зависимости μ от r")
    canvas.draw()

def toggle_show_points():
    global show_points
    show_points = var_show_points.get()
    update_plot(auto=False)

def create_table(num_points):
    global table
    for child in table_frame.winfo_children():
        child.destroy()
    columns = ("r", "μ")
    table = ttk.Treeview(table_frame, columns=columns, show="headings")
    table.heading("r", text="r")
    table.heading("μ", text="μ")
    table.grid(row=0, column=0, sticky=tk.NSEW)
    for i in range(1, num_points + 1):
        table.insert("", "end", values=(i, 0.0))
    table.bind("<Double-1>", on_double_click)
    scrollbar = ttk.Scrollbar(table_frame, orient="vertical", command=table.yview)
    table.configure(yscrollcommand=scrollbar.set)
    scrollbar.grid(row=0, column=1, sticky=tk.NS)

def on_double_click(event):
    region = table.identify_region(event.x, event.y)
    if region == "cell":
        column = table.identify_column(event.x)
        row = table.identify_row(event.y)
        if row and column:
            cell_value = table.item(row, "values")[int(column[1:]) - 1]
            entry_edit = ttk.Entry(tab1)
            entry_edit.insert(0, cell_value)
            x, y, width, height = table.bbox(row, column)
            entry_edit.place(x=x + 55, y=y + 470, width=width, height=height)

            def save_edit(event):
                new_value = entry_edit.get()
                values = list(table.item(row, "values"))
                values[int(column[1:]) - 1] = new_value
                table.item(row, values=values)
                entry_edit.destroy()

            entry_edit.bind("<Return>", save_edit)
            entry_edit.bind("<FocusOut>", lambda e: entry_edit.destroy())

def load_from_file():
    file_path = filedialog.askopenfilename(filetypes=[("Text files", "*.txt"), ("All files", "*.*")])
    if file_path:
        with open(file_path, "r") as file:
            lines = file.readlines()
            for i, line in enumerate(lines):
                if i < len(table.get_children()):
                    values = line.strip().split()
                    if len(values) == 2:
                        try:
                            r = float(values[0])
                            mu = float(values[1])
                            table.item(table.get_children()[i], values=[r, mu])
                        except ValueError:
                            messagebox.showerror("Ошибка", f"Некорректные данные в строке {i + 1}: {line.strip()}")
                            return
                else:
                    break

# Создание интерфейса
root = tk.Tk()
root.title("Аппроксиматор")
root.geometry("1200x1000")
var_mu = tk.IntVar(value=1)

tab_control = ttk.Notebook(root)
tab1 = ttk.Frame(tab_control)


# Вкладка 1
label_title_tab1 = tk.Label(tab1, text="Ввод исходных данных", font=("Arial", 14))
label_title_tab1.grid(row=0, column=0, columnspan=2, pady=10)

fields_tab1 = [
    ("Модуль сдвига однородного материала μ₀", "mu0"),
    ("Граница а", "a"),
    ("Граница b", "b")
]

entries = {}
for i, (field, name) in enumerate(fields_tab1):
    label = tk.Label(tab1, text=field)
    label.grid(row=i + 1, column=0, sticky=tk.W, padx=5, pady=5)
    entry = tk.Entry(tab1)
    entry.grid(row=i + 1, column=1, padx=5, pady=5)
    entries[name] = entry

entry_a = entries['a']
entry_b = entries['b']

tab2 = Tab2(tab_control, data_mu1, data_mu2, {
    'var_mu': var_mu,
    'mu0': entries['mu0'].get()
})
tab_control.add(tab1, text="Вкладка 1")
tab_control.add(tab2.tab, text="Вкладка 2")
tab_control.pack(expand=1, fill="both")

frame_mu = tk.Frame(tab1)
frame_mu.grid(row=len(fields_tab1) + 1, column=0, columnspan=2, pady=10)
tk.Radiobutton(frame_mu, text="μ₁", variable=var_mu, value=1, command=switch_mu).pack(side=tk.LEFT)
tk.Radiobutton(frame_mu, text="μ₂", variable=var_mu, value=2, command=switch_mu).pack(side=tk.LEFT)

combo_mu1 = ttk.Combobox(tab1, values=[
    "μ(r) = μ₀(1 - k _exp(-α_ r))",
    "μ(r) = μ₀(1 + k_ ((b - r) / (b - a))^n)",
    "μ(r) = μ₀[1 + (k-1)(a/r)^n]",
    "μ(r) = m(1 + k e^{-ar} (cos(βr + ω) + sin(βr + ω)))"
])
combo_mu1.grid(row=len(fields_tab1) + 2, column=0, columnspan=2, pady=5)
combo_mu1.current(0)

combo_mu2 = ttk.Combobox(tab1, values=[
    "μ(r) = μ₀(1 - k exp(-α r))",
    "μ(r) = μ₀(1 + k_ ((b - r) / (b - a))^n)",
    "μ(r) = μ₀[1 + (k-1)(a/r)^n]",
    "μ(r) = m(1 + k e^{-ar} (cos(βr + ω) + sin(βr + ω)))"
])
combo_mu2.grid(row=len(fields_tab1) + 3, column=0, columnspan=2, pady=5)
combo_mu2.current(0)
combo_mu2.grid_forget()

label_approx = tk.Label(tab1, text="Узловые точки для аппроксимации модуля сдвига")
label_approx.grid(row=len(fields_tab1) + 4, column=0, columnspan=2, pady=10)

num_points_var = tk.IntVar(value=10)
num_points_combobox = ttk.Combobox(tab1, textvariable=num_points_var, values=list(range(5, 21, 5)))
num_points_combobox.grid(row=len(fields_tab1) + 5, column=0, columnspan=2, pady=5)
num_points_combobox.bind("<<ComboboxSelected>>", lambda event: create_table(num_points_var.get()))

table_frame = tk.Frame(tab1)
table_frame.grid(row=len(fields_tab1) + 6, column=0, columnspan=2, pady=10)

fig, ax = plt.subplots()
canvas = FigureCanvasTkAgg(fig, master=tab1)
canvas.get_tk_widget().grid(row=0, column=2, rowspan=10, padx=10, pady=10, sticky=tk.NSEW)

# После создания canvas:
canvas.mpl_connect("scroll_event", on_mousewheel)
canvas.mpl_connect("button_press_event", on_press)
canvas.mpl_connect("button_release_event", on_release)
canvas.mpl_connect("motion_notify_event", on_motion)
# toolbar_frame = tk.Frame(tab1)
# toolbar_frame.grid(row=9, column=2, padx=10, pady=2, sticky=tk.NSEW)
# toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
# toolbar.update()

button_update = tk.Button(tab1, text="Обновить график", command=lambda: [save_current_mu_data(), update_plot(auto=False)])
button_update.grid(row=len(fields_tab1) + 7, column=0, columnspan=2, pady=10)

button_auto_update = tk.Button(tab1, text="Автоматическое отображение графика", command=lambda: [save_current_mu_data(), update_plot(auto=True)])
button_auto_update.grid(row=len(fields_tab1) + 8, column=0, columnspan=2, pady=10)

button_save = tk.Button(tab1, text="Сохранить", command=save_data)
button_save.grid(row=len(fields_tab1) + 9, column=0, columnspan=2, pady=10)

button_load = tk.Button(tab1, text="Считать из файла", command=load_from_file)
button_load.grid(row=len(fields_tab1) + 10, column=0, columnspan=2, pady=10)

button_clear_plot = tk.Button(tab1, text="Очистить график", command=clear_plot)
button_clear_plot.grid(row=len(fields_tab1) + 11, column=0, columnspan=2, pady=10)

# Поля для настройки осей
label_xmin = tk.Label(tab1, text="X min:")
label_xmin.grid(row=len(fields_tab1) + 12, column=0, sticky=tk.W, padx=5, pady=5)
entry_xmin = tk.Entry(tab1)
entry_xmin.grid(row=len(fields_tab1) + 12, column=1, padx=5, pady=5)
entry_xmin.insert(0, "0")

label_xmax = tk.Label(tab1, text="X max:")
label_xmax.grid(row=len(fields_tab1) + 13, column=0, sticky=tk.W, padx=5, pady=5)
entry_xmax = tk.Entry(tab1)
entry_xmax.grid(row=len(fields_tab1) + 13, column=1, padx=5, pady=5)
entry_xmax.insert(0, "3")

label_ymin = tk.Label(tab1, text="Y min:")
label_ymin.grid(row=len(fields_tab1) + 14, column=0, sticky=tk.W, padx=5, pady=5)
entry_ymin = tk.Entry(tab1)
entry_ymin.grid(row=len(fields_tab1) + 14, column=1, padx=5, pady=5)
entry_ymin.insert(0, "0")

label_ymax = tk.Label(tab1, text="Y max:")
label_ymax.grid(row=len(fields_tab1) + 15, column=0, sticky=tk.W, padx=5, pady=5)
entry_ymax = tk.Entry(tab1)
entry_ymax.grid(row=len(fields_tab1) + 15, column=1, padx=5, pady=5)
entry_ymax.insert(0, "1")

label_xstep = tk.Label(tab1, text="X step:")
label_xstep.grid(row=len(fields_tab1) + 16, column=0, sticky=tk.W, padx=5, pady=5)
entry_xstep = tk.Entry(tab1)
entry_xstep.grid(row=len(fields_tab1) + 16, column=1, padx=5, pady=5)
entry_xstep.insert(0, "0.333")

label_ystep = tk.Label(tab1, text="Y step:")
label_ystep.grid(row=len(fields_tab1) + 17, column=0, sticky=tk.W, padx=5, pady=5)
entry_ystep = tk.Entry(tab1)
entry_ystep.grid(row=len(fields_tab1) + 17, column=1, padx=5, pady=5)
entry_ystep.insert(0, "0.25")

# Добавление чекбокса "отображение точек"
var_show_points = tk.BooleanVar(value=True)
checkbutton_points = tk.Checkbutton(tab1, text="Отображение точек", variable=var_show_points, command=toggle_show_points)
checkbutton_points.grid(row=len(fields_tab1) + 18, column=0, columnspan=2, pady=10)

create_table(num_points_var.get())
root.mainloop()