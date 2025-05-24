from sqlalchemy import create_engine
import pandas as pd
from sklearn 

sql_engine = create_engine("postgresql://chrismitsacopoulos:password@localhost/pasteurdb")
df = pd.read_sql_query("SELECT * FROM training_dataset", sql_engine)

#Using regression, it resembles any simple function; designate X and Y (desired output parameters).
X_training_set = df["h3_chain", "l3_chain", "h3_gravy", "l3_gravy", "h3_pi", "l3_pi", "h3_n_gylcosylation_sites", "h3_o_gylcosylation_sites", "l3_n_gylcosylation_sites", "l3_o_gylcosylation_sites", ]
print(len(df.columns))