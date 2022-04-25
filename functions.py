import scipy as sp
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

df = pd.read_csv('tables/antconst.csv')
df_prop = pd.read_csv('tables/thermophysprop.csv')

names_antoine = list(df['Name'])
names_physprop = list(df_prop['Name'])

vap_press = 0

def calculate_antoine(species, temp):
    # print('Species are :', names_antoine)
    index = df[df['Name'] == species].index.values
    ABC = df.iloc[index, 2:5].values.tolist()
    try:
        A = float(ABC[0][0])
        B = float(ABC[0][1])
        C = float(ABC[0][2])

    except ValueError:
        print('invalid input')
        # answer.configure(text= 'invalid input')
    else:
        vap_press = np.exp(A - (B / (temp + C)))
        return vap_press
        # answer.configure(text=f'{round(vap_press, 3)} kPa')

def get_species_properties():
    print(names_physprop)
    species_name = input('type species name')
    index_2 = df_prop[df_prop['Name'] == species_name].index.values
    species_values = df_prop.iloc[index_2, 1:].values.tolist()
    print(['Molar mass','ω','Tc/K','Pc/bar','Zc','Vc·cm3·mol−1','Tn/K'])
    print(species_values)

def calculate_lee_kessler(species_name,temperature,pressure):
    #species_name = 'Water' #input('species name')
    index_2 = df_prop[df_prop['Name'] == species_name].index.values
    species_values = df_prop.iloc[index_2, 2:].values.tolist()
    #pressure = #float(input('type in pressure'))
    #temperature = #float(input('type in temperature'))
    temp_r = temperature / species_values[0][3]
    press_r = pressure / species_values[0][4]
    if temp_r <= 1:
        df_kess_o = pd.read_csv('tables/Kessler_HRo.csv')
        df_kess_1 = pd.read_csv('tables/Kessler_HR1.csv')
        value_o = df_kess_o[press_r, temp_r]
        value_1 = df_kess_o[press_r, temp_r]
    elif temp_r >=1:
        df_kess_o = pd.read_csv('tables/Kessler_HRo_cont.csv')
        df_kess_1 = pd.read_csv('tables/Kessler_HR1_cont.csv')
        value_o = df_kess_o[press_r, temp_r]
        value_1 = df_kess_o[press_r, temp_r]
    else:
        print('error')
    return value_o, value_1



def calculate_bubbl(species_name_1,species_name_2,x_1):
    # Returns dew point and bubble point at a specific molar fraction for the solution
    species_name_1 = "Water"
    species_name_2 = "Ethanol"
    index_spec_1 = df_prop[df_prop['Name'] == species_name_1].index.values
    index_spec_2 = df_prop[df_prop['Name'] == species_name_2].index.values
    species_values_1 = df_prop.iloc[index_spec_1, 2:5].values.tolist()
    species_values_2 = df_prop.iloc[index_spec_2, 2:5].values.tolist()
    p_vap_spec_1 = np.exp(species_values_1[0] - (species_values_1[1] / (temp + species_values_1[2])))
    p_vap_spec_2 = np.exp(species_values_2[0] - (species_values_2[1] / (temp + species_values_2[2])))

    #bubbl_temp =

def calculate_dew(species_name_1, species_name_2):
    # Returns dew point and bubble point at a specific molar fraction for the solution
    species_name_1 = "Water"
    species_name_2 = "Ethanol"
    index_spec_1 = df_prop[df_prop['Name'] == species_name_1].index.values
    index_spec_2 = df_prop[df_prop['Name'] == species_name_2].index.values
    species_values_1 = df_prop.iloc[index_spec_1, 2:5].values.tolist()
    species_values_2 = df_prop.iloc[index_spec_2, 2:5].values.tolist()
    p_vap_spec_1 = np.exp(species_values_1[0] - (species_values_1[1] / (temp + species_values_1[2])))
    p_vap_spec_2 = np.exp(species_values_2[0] - (species_values_2[1] / (temp + species_values_2[2])))

    dew_point = 1 / ((y_1 / p_vap_spec_1) + ((1 - y_1) / p_vap_spec_2))


def q_line(x,z,q):
    y = (q/(q-1))*x - z/(q-1)
    return y

def rectifying_line(R,x,xd):
    y = (R/(R+1))*x+((1/(R+1))*xd)
    return y 

