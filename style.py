import tkinter as tk
from tkinter import ttk
import customtkinter as ctk


def configure_styles(root):
    """Настройка стилей для ttk и tk виджетов"""
    # Настройка ttk.Style
    style = ttk.Style()
    style.theme_use('clam')  # Используем тему 'clam'

    # Стиль для кнопок ttk
    style.configure('TButton',
                    background='lightblue',
                    foreground='black',
                    font=('Arial', 12),
                    padding=5,
                    borderwidth=2,
                    relief="ridge")

    # Стиль для полей ввода ttk
    style.configure('TEntry',
                    fieldbackground='lightyellow',
                    foreground='darkblue',
                    font=('Arial', 12),
                    padding=3)

    # Глобальные настройки для tk.Entry через опции
    root.option_add('*Entry*Background', 'lightyellow')
    root.option_add('*Entry*Foreground', 'darkblue')
    root.option_add('*Entry*Font', ('Arial', 12))

    # Настройка tk.Button
    tk.Button.configure(root,
                        bg='lightgreen',
                        fg='black',
                        font=('Arial', 10, 'bold'),
                        activebackground='darkgreen',
                        activeforeground='white',
                        borderwidth=2,
                        relief="groove")

    # Настройка меток
    root.option_add('*Label*Font', ('Arial', 12))
    root.option_add('*Label*Background', 'white')

    # Настройка таблиц (Treeview)
    style.configure('Treeview',
                    background='lightgray',
                    fieldbackground='lightgray',
                    foreground='black',
                    font=('Arial', 10))
    style.map('Treeview', background=[('selected', 'blue')])

def configure_borderRad(root):
    ctk.set_appearance_mode("System")
    ctk.set_default_color_theme("blue")

def apply_custom_radius(widget, radius=10):
    """Функция для эмуляции border-radius (требуется Pillow)"""
    # Для реализации border-radius потребуется Pillow и customtkinter
    pass  # Реализуйте через customtkinter или другие методы