import pandas as pd

df = pd.read_csv('result/result.csv')
columns = df.columns
return_columns = []
for column in columns:
    if (column[:2] == 'As') or (column[:1] == 'e') or (column in ['Tmax', 'Tmin']):
        return_columns.append(column)
        
df_display = df[return_columns]