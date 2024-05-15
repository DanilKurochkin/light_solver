import numpy as np
import yaml
from classes import light_solver as ls
from classes.cyclogram import Cyclogram
from constans import labels, args

with open("inputs/light_solve.yaml") as f:
    data = yaml.full_load(f)

Tmin = data['Tmin'] # 0 *c
Tmax = data['Tmax'] # +30 *c
width = data['width']
Lx = data['Lx']
Ly = data['Ly']
Lz = data['Lz']
c = data['c']
p = data['p']
R = data['R']
orbit_height = data['orbit_height']
As = data['As']
e = data['e']

cyclograms = []
for item in data['generation_heat']:
    heat = np.array(item)[0]
    time = np.array(item)[1]
    cyclograms.append(Cyclogram(heat, time))

new_columns_for_T = []
for i in range(len(cyclograms)):
    for arg in args:
        for label in labels:
            new_columns_for_T.append(f"{label}_{arg}_{i+1}")
    
    for j in range(cyclograms[i].peaks_amount):
        for lable in labels:
            new_columns_for_T.append(f"{label}_peak_{j+1}_{i+1}")

solver = ls.Solver(
    width=width, Lx=Lx, Ly=Ly, Lz=Lz,
    c=c, p=p, R=R,
    orbit_height=orbit_height,
    cyclograms=cyclograms)

print(solver.Qav_As)
print(solver.Qav_e)

x=solver.temperature_diapozone(
    As=[0.1,
        0.8,
        0.9,
        0.1,
        0.9,
        0.9],
    e= [0.1,
        0.1,
        0.3,
        0.5,
        0.6,
        0.7],
    n=1
)

print(x[0:6])
print(x[12:18])
# df, As_cols, e_cols = solver.construct_df(
#     As, e, Tmin, Tmax
# )
# print(f"Finded: {len(df)} combination")

# n = len(df)
# df[new_columns_for_T] = np.vectorize(
#     solver.temperature_diapozone,
#     signature='(m), (m), ()->(n)')(df[As_cols], df[e_cols], n)

# df['Tmin'] = df[new_columns_for_T].min(axis=1)
# df['Tmax'] = df[new_columns_for_T].max(axis=1)
# df.to_csv('result/result.csv')