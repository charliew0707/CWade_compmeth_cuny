# Computational Methods in Physics

Computational Methods in Physics [Class # 54168] – PHYS 78100
Graduate Center Campus Tu, Th, 9:45 am – 11:35 pm Room #7395 4.000 Credits Ariyeh Maller

Course repository for Phys 78100 at CUNY, Fall 2026. Taught by Professor Maller, following *Computational Physics* by Mark Newman.

The course is an introduction to numerical techniques used in physics and astrophysics — topics include numerical integration and differentiation, solving differential equations, Fourier transforms, Monte Carlo methods, and high performance computing, all implemented in Python.

---

## What's in here

**`lecture_notes/`** — Jupyter notebooks (`.ipynb`) from each class session. These are interactive notebooks that mix code, plots, and written explanations, worked through during lecture.

**`exercises/`** — In-class coding exercises done during lecture, also as Jupyter notebooks.

**`homework/`** — Homework assignments written in VSCode and submitted via Git. Python scripts (`.py`) or notebooks depending on the assignment.

**`final_project/`** — Final project code and presentation materials. The project involves choosing a physics problem and solving it numerically using techniques from the course.

---

## Setup

Dependencies are managed with [uv](https://github.com/astral-sh/uv). To get started:

```bash
uv sync
```

Or install packages individually:
```bash
uv add numpy scipy matplotlib sympy jupyterlab
```

Launch Jupyter for lecture notes and exercises:
```bash
uv run jupyter lab
```

---

*Charles Wade | charliew0707@gmail.com*
