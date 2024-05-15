import pandas as pd
import streamlit as st
from constans import labels, args
from classes import reader

df_display = reader.df_display

st.title("Выбор покрытия панели")

Tmin = st.number_input(
    label="Минимальная температура", value=None
)
Tmax = st.number_input(
    label="Максимальная температура", value=None
)

As_col, e_col = st.columns(2)

As = []
e = []
for i in range(6):
    As.append(As_col.number_input(
        label=f"As для {labels[i]}",
        min_value=0., max_value=1., key=f'As_{i}', value=None
        )
    )
    
    e.append(e_col.number_input(
        label=f"e для {labels[i]}",
        min_value=0., max_value=1., key=f'e_{i}', value=None
        )
    )

if st.button("Подобрать подходящие"):
    df_temp = df_display[(df_display["Tmax"] <= Tmax) & (df_display["Tmin"] >= Tmin)]
    for i in range(6):
        if As[i] is not None:
            df_temp = df_temp[df_temp[f'As_{i}'] == As[i]]
        if e[i] is not None:
            df_temp = df_temp[df_temp[f'e_{i}'] == e[i]]
    
    st.write(f"Найдено {len(df_temp)} комбинаций")
    st.dataframe(df_temp)