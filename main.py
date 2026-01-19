import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import sympy as sp
from scipy.stats import binom, poisson
import re


class UltraCalcV23:
    def __init__(self, root):
        self.root = root
        self.root.title("Professional Math & Statistics Suite v2.3")
        self.root.geometry("1300x950")
        self.root.configure(bg="#1e1e2e")

        # Grid weight configuration for full-screen scaling
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)

        # Global Font Scaling (+6 increase)
        self.font_main = ("Segoe UI", 16)
        self.font_header = ("Segoe UI", 18, "bold")
        self.font_mono = ("Consolas", 16)

        self.style = ttk.Style()
        self.style.theme_use('clam')
        self.style.configure("TNotebook", background="#1e1e2e", borderwidth=0)
        self.style.configure("TNotebook.Tab", font=("Segoe UI", 14), padding=[10, 5])
        self.style.configure("TFrame", background="#f8f9fa")
        self.style.configure("TLabel", background="#f8f9fa", font=self.font_main)
        self.style.configure("TLabelframe", background="#f8f9fa")
        self.style.configure("TLabelframe.Label", font=self.font_header, background="#f8f9fa")

        self.main_container = tk.Frame(self.root, bg="#1e1e2e")
        self.main_container.grid(row=0, column=0, sticky="nsew")
        self.main_container.columnconfigure(1, weight=1)
        self.main_container.rowconfigure(0, weight=1)

        self.setup_symbol_sidebar()

        self.notebook = ttk.Notebook(self.main_container)
        self.notebook.grid(row=0, column=1, sticky="nsew", padx=15, pady=15)

        self.setup_algebra_tab()
        self.setup_statistics_tab()

        # Bottom Export Bar
        export_frame = tk.Frame(self.root, bg="#1e1e2e")
        export_frame.grid(row=1, column=0, sticky="ew")
        tk.Button(export_frame, text="📥 EXPORT DATA",
                  command=self.export_results, bg="#4CAF50", fg="white",
                  font=("Arial", 14, "bold"), pady=12).pack(fill="x", padx=25, pady=15)

    def setup_symbol_sidebar(self):
        sidebar = tk.Frame(self.main_container, bg="#2d2d3d", width=150)
        sidebar.grid(row=0, column=0, sticky="ns", padx=5, pady=10)
        tk.Label(sidebar, text="SYMBOLS", bg="#2d2d3d", fg="white", font=("Arial", 12, "bold")).pack(pady=15)

        symbols = [("xⁿ", "^"), ("√", "sqrt()"), ("π", "pi"), ("!", "!"), ("nCr", "C"), ("i", "I"), ("(", "("),
                   (")", ")")]
        for text, code in symbols:
            tk.Button(sidebar, text=text, command=lambda c=code: self.insert_symbol(c),
                      bg="#3d3d5d", fg="white", relief="flat", width=10, font=("Arial", 12), pady=8).pack(pady=4,
                                                                                                          padx=12)

    def insert_symbol(self, char):
        focused = self.root.focus_get()
        if isinstance(focused, (tk.Entry, tk.Text)):
            focused.insert(tk.INSERT, char)

    def preprocess_input(self, text):
        text = text.replace('^', '**')
        # Fix implicit multiplication: 2x -> 2*x, (x)2 -> (x)*2, (x)(x) -> (x)*(x)
        text = re.sub(r'(\d)([a-zA-Z\(])', r'\1*\2', text)
        text = re.sub(r'(\))(\d|[a-zA-Z\(])', r'\1*\2', text)
        return text

    def setup_algebra_tab(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="   Algebra & Polynomials   ")
        tab.columnconfigure(0, weight=1)
        tab.rowconfigure(2, weight=1)

        # Solver Section
        p_frame = ttk.LabelFrame(tab, text=" Polynomial Solver & Grapher ")
        p_frame.grid(row=0, column=0, sticky="ew", padx=25, pady=15)
        p_frame.columnconfigure(0, weight=1)

        self.poly_input = ttk.Entry(p_frame, font=self.font_mono)
        self.poly_input.grid(row=0, column=0, sticky="ew", padx=15, pady=20)
        self.poly_input.insert(0, "(2x^4)+(2x)-1")

        btn_row = ttk.Frame(p_frame)
        btn_row.grid(row=0, column=1, padx=10)
        ttk.Button(btn_row, text="SOLVE ROOTS", command=self.solve_poly).pack(side="left", padx=10)
        ttk.Button(btn_row, text="GRAPH", command=self.graph_poly).pack(side="left", padx=10)

        # Binomial Section
        b_frame = ttk.LabelFrame(tab, text=" Binomial Expansion ")
        b_frame.grid(row=1, column=0, sticky="ew", padx=25, pady=15)
        row = ttk.Frame(b_frame);
        row.pack(pady=15)
        self.a_in = ttk.Entry(row, width=10, font=self.font_mono);
        self.a_in.insert(0, "x")
        self.b_in = ttk.Entry(row, width=10, font=self.font_mono);
        self.b_in.insert(0, "y")
        self.n_in = ttk.Entry(row, width=5, font=self.font_mono);
        self.n_in.insert(0, "3")

        ttk.Label(row, text=" ( ", font=("Arial", 20)).pack(side="left")
        self.a_in.pack(side="left", padx=5)
        ttk.Label(row, text=" + ", font=("Arial", 20)).pack(side="left")
        self.b_in.pack(side="left", padx=5)
        ttk.Label(row, text=" ) ^ ", font=("Arial", 20)).pack(side="left")
        self.n_in.pack(side="left", padx=5)
        ttk.Button(b_frame, text="EXPAND", command=self.expand_bin).pack(pady=10)

        # Output Text (Scalable)
        self.alg_txt = tk.Text(tab, font=self.font_mono, bg="#ffffff", padx=15, pady=15, undo=True)
        self.alg_txt.grid(row=2, column=0, sticky="nsew", padx=25, pady=20)

    def solve_poly(self):
        try:
            expr_str = self.preprocess_input(self.poly_input.get())
            x = sp.Symbol('x')
            expr = sp.sympify(expr_str)
            roots = sp.solve(expr, x)

            self.alg_txt.delete('1.0', tk.END)
            self.alg_txt.insert(tk.END, f"EQUATION: {sp.pretty(expr)}\n" + "─" * 70 + "\n")
            self.alg_txt.insert(tk.END, f"Total Roots Found: {len(roots)}\n\n")

            for i, r in enumerate(roots):
                # evalf() ensures complex/imaginary parts are displayed numerically
                self.alg_txt.insert(tk.END, f"Root x{i + 1}: {r.evalf()}\n")
        except Exception as e:
            messagebox.showerror("Math Error", f"Check syntax: {e}")

    def graph_poly(self):
        try:
            expr_str = self.preprocess_input(self.poly_input.get())
            x_sym = sp.Symbol('x')
            expr = sp.sympify(expr_str)
            f = sp.lambdify(x_sym, expr, 'numpy')

            graph_win = tk.Toplevel(self.root)
            graph_win.title("Interactive Function Graph")
            graph_win.geometry("900x700")

            fig, ax = plt.subplots(figsize=(8, 6), dpi=100)
            x_range = np.linspace(-10, 10, 1000)
            y_range = f(x_range)

            ax.plot(x_range, y_range, color='#e74c3c', linewidth=2.5, label=f"y = {expr}")
            ax.axhline(0, color='black', linewidth=1.2)
            ax.axvline(0, color='black', linewidth=1.2)
            ax.grid(True, which='both', linestyle='--', alpha=0.6)
            ax.set_xlabel("x", fontsize=14)
            ax.set_ylabel("f(x)", fontsize=14)
            ax.legend(fontsize=12)

            canvas = FigureCanvasTkAgg(fig, master=graph_win)
            canvas.get_tk_widget().pack(fill="both", expand=True)

            # Toolbar for zooming and panning
            toolbar = NavigationToolbar2Tk(canvas, graph_win)
            toolbar.update()
            canvas.draw()
        except:
            messagebox.showerror("Graphing Error", "Check the polynomial expression.")

    def setup_statistics_tab(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="   Statistical Distributions   ")
        tab.columnconfigure(1, weight=3)
        tab.rowconfigure(0, weight=1)

        ctrl = ttk.Frame(tab)
        ctrl.grid(row=0, column=0, sticky="nsew", padx=20, pady=20)

        self.st_model = ttk.Combobox(ctrl, values=["Binomial", "Poisson"], state="readonly", font=self.font_main)
        self.st_model.pack(fill="x", pady=15)
        self.st_model.bind("<<ComboboxSelected>>", self.draw_stat_inputs)

        self.st_inputs = ttk.Frame(ctrl);
        self.st_inputs.pack(fill="x", pady=15)
        self.st_res = tk.StringVar(value="P(X) = 0.0000")
        tk.Label(ctrl, textvariable=self.st_res, font=("Arial", 18, "bold"), fg="#d35400").pack(pady=30)
        ttk.Button(ctrl, text="CALCULATE", command=self.run_stats).pack(fill="x", pady=10)

        self.fig_s, self.ax_s = plt.subplots(figsize=(7, 5))
        self.canvas_s = FigureCanvasTkAgg(self.fig_s, master=tab)
        self.canvas_s.get_tk_widget().grid(row=0, column=1, sticky="nsew", padx=20, pady=20)

    def draw_stat_inputs(self, e=None):
        for w in self.st_inputs.winfo_children(): w.destroy()

        # Increased font for dropdown and entries
        ttk.Label(self.st_inputs, text="Inequality:").grid(row=0, column=0, sticky="w", pady=8)
        self.st_sym = ttk.Combobox(self.st_inputs, values=["x ≤ k", "x ≥ k", "x = k"], font=self.font_main, width=12)
        self.st_sym.grid(row=0, column=1, pady=8, padx=10)

        ttk.Label(self.st_inputs, text="Value k:").grid(row=1, column=0, sticky="w", pady=8)
        self.st_k = ttk.Entry(self.st_inputs, font=self.font_main, width=14);
        self.st_k.grid(row=1, column=1, pady=8, padx=10)

        if self.st_model.get() == "Binomial":
            ttk.Label(self.st_inputs, text="Trials (n):").grid(row=2, column=0, sticky="w", pady=8)
            self.st_n = ttk.Entry(self.st_inputs, font=self.font_main, width=14);
            self.st_n.grid(row=2, column=1, pady=8, padx=10)
            ttk.Label(self.st_inputs, text="Prob (p):").grid(row=3, column=0, sticky="w", pady=8)
            self.st_p = ttk.Entry(self.st_inputs, font=self.font_main, width=14);
            self.st_p.grid(row=3, column=1, pady=8, padx=10)
        else:
            ttk.Label(self.st_inputs, text="Rate (λ):").grid(row=2, column=0, sticky="w", pady=8)
            self.st_lam = ttk.Entry(self.st_inputs, font=self.font_main, width=14);
            self.st_lam.grid(row=2, column=1, pady=8, padx=10)

    def expand_bin(self):
        try:
            a = sp.sympify(self.preprocess_input(self.a_in.get()))
            b = sp.sympify(self.preprocess_input(self.b_in.get()))
            n = int(self.n_in.get())
            res = sp.expand((a + b) ** n)
            self.alg_txt.delete('1.0', tk.END)
            self.alg_txt.insert(tk.END, f"BINOMIAL EXPANSION of ({a} + {b})^{n}:\n\n")
            self.alg_txt.insert(tk.END, sp.pretty(res, use_unicode=True))
        except:
            pass

    def run_stats(self):
        try:
            self.ax_s.clear()
            model, sym, k = self.st_model.get(), self.st_sym.get(), int(self.st_k.get())
            if model == "Binomial":
                n, p = int(self.st_n.get()), float(self.st_p.get())
                dist = binom(n, p);
                x_vals = np.arange(0, n + 1)
            else:
                lam = float(self.st_lam.get());
                dist = poisson(lam)
                x_vals = np.arange(0, int(lam * 3) + 5)

            prob = dist.cdf(k) if "≤" in sym else (1 - dist.cdf(k - 1) if "≥" in sym else dist.pmf(k))
            self.st_res.set(f"P({sym}) = {prob:.6f}")
            self.ax_s.bar(x_vals, dist.pmf(x_vals), color='#3498db', alpha=0.8, edgecolor='black')
            self.ax_s.set_title(f"{model} Distribution", fontsize=16)
            self.canvas_s.draw()
        except:
            pass

    def export_results(self):
        content = self.alg_txt.get("1.0", tk.END)
        fpath = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text File", "*.txt")])
        if fpath:
            with open(fpath, "w", encoding="utf-8") as f: f.write(content)
            messagebox.showinfo("Success", "All results exported successfully.")


if __name__ == "__main__":
    root = tk.Tk()
    app = UltraCalcV23(root)
    root.mainloop()