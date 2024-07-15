import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#Constants 
fsy = 500     # [MPa] Steel Yield Strength
Es = 200000         # [MPa] Steel Modulus (Constant)
ec = 0.003          # [MPa] Concrete Maximum strain (MPa)


def calculate_column_design(col_sect, col_type, fc, b, h, dbar, xbar, ybar, l, ke2, ke3, cover, fy, tie_db, tie_spc):
    def Min_TieSize(db):
        if db <= 20:
            Min_TieSize = 6
        elif db <= 24 or db <= 28:
            Min_TieSize = 10
        elif db <= 32 or db <= 36:
            Min_TieSize = 12
        else:  
            Min_TieSize = 16
        return Min_TieSize

    tie = max(Min_TieSize(dbar), tie_db)        # Ensure Tie Size complies with Code

    Max_TieSpacing = min(b, h, 15*dbar)
    if fc > 50: 
        Max_TieSpacing = min(0.8*min(b, h), min(b,h), 300)

    tie_spc = min(tie_spc, Max_TieSpacing)
    if tie_spc <= 0: 
        tie_spc = Max_TieSpacing

    def col_Areas(sect, b, D):     # [mm^4] Second Moment of inertia
        match sect:
            case 'RECT':
                Ag = b * D                          # [mm²] Column Area
                Ig = b * D ** 3 / 12                # [mm^4] Second Moment of inertia
            case 'CIRC':
                Ag = math.pi * (D/2) **2
                Ig = math.pi /4 * (D/2) ** 4
        rgy = (Ig/ Ag) ** 0.5                       # [mm] Radius of Gyration
        return Ag, Ig, rgy

    Col3 = col_Areas(col_sect, b, h)
    Ag = Col3[0]
    I3 = Col3[1]
    r3 = Col3[2]

    Col2 = col_Areas(col_sect, h, b)
    I2 = Col2[1]
    r2 = Col2[2]

    M3min = h * 0.05                                 # [kNm] Minimum Moment (y-dirn)
    M2min = b * 0.05                                 # [kNm] Minimum Moment (x-dirn)

    deff = h - cover - tie - dbar/2     # [mm] Distance from uppermost steel centre to bottom edge 

    dox = b - cover - tie - dbar/2      # [mm] Distance from uppermost steel centre to bottom edge (X-Dirn) 
    doy = h - cover - tie - dbar/2      # [mm] Distance from uppermost steel centre to bottom edge (y-Dirn)

    dprime = cover + tie + dbar/2       # [mm] Distance from top of concrete to centre of reinforcement

    if fc <= 28:
        beta1 = 0.85
    else:
        beta1 = max(0.85 - 0.05 * (fc - 28) / 7, 0.65)

    alpha1 = min(0.85, max(1 - 0.003 * fc, 0.72))

    def alpha2(fc, section):
        alpha2 = max(0.85 - 0.0015 * fc, 0.67)
        
        if section == 'CIRC':
            alpha2 = alpha2 * 0.95

        return alpha2

    gamma = max(0.97-0.0025 * fc, 0.67)

    match col_sect:
        case 'RECT':
            Bar_total = 2 * (xbar + ybar-2)
        case 'CIRC':
            Bar_total = ybar

    Ast = (math.pi *dbar**2/4) * Bar_total      # [mm²] Steel area at tension side
    Asc = Ast                                   # [mm²] Steel area at compression side

    ey = fy / Es        # Steel strain at yield point

    def stress_strain(c):
        es = ec * (deff -c)/c       # Steel strain at neutral axis
        fs = es * Es
        esprime = ec *(c-dprime)/c  # Steel strain at compression side
        fsprime = esprime * Es      #tensile stress at compression side 
        a = beta1* c                # Compressive stress block depth
        return a, es, fs, esprime, fsprime

    def strength_factor():
        phi_factor = 0.65 + 0.25 * (es - ey)/(0.005-ey)
        if phi_factor >= 0.9:
            phi, classify = 0.9, 'Tension-controlled'
        elif phi_factor <= 0.65:
            phi, classify = 0.65, 'Compression-controlled'
        else:
            phi, classify = phi_factor, 'Transition'
        return phi, classify

    def forces_moments():
        T = Ast * min(fy,fs) /1000                              # [kN] Tension force by steel at tension side
        Cc = 0.85* fc * a * b /1000                             # [kN] Compressive force by concrete
        Cs = Asc * (min(fy, fsprime) - (0.85*fc)) /1000         # [kN] Compressive force by steel at compression side

        Pn = Cc + Cs - T                                        # [kN] Nominal Axial load capacity 
        Pu = phi * Pn                                           # [kN] Ultimate axial load capacity

        ecc = (Cc*(deff - a/2) + Cs*(deff - dprime)) / Pn - (h/2 - dprime)
        Mn = Pn * ecc /1000                             # Nominal Moment Capacity
        Mu = phi * Mn                                   # Ultimate Moment Capacity

        return ecc, Pn, Pu, Mn, Mu

    nom_load = []
    ult_load = []
    nom_moment = []
    ult_moment = []
    eccentricity = []
    phi_factor = []

    c_values = np.arange(1000, 50, -5)
    for c in c_values:
        a, es, fs, esprime, fsprime = stress_strain(c)
        phi, classify = strength_factor()
        ecc, Pn, Pu, Mn, Mu = forces_moments()
        if ecc <= 1.5*h:
            nom_load.append(round(Pn))
            ult_load.append(round(Pu))
            nom_moment.append(round(Mn))
            ult_moment.append(round(Mu))
            eccentricity.append(round(ecc))
            phi_factor.append(round(phi,2))

    dict = {'ecc(mm)': eccentricity, 'Pn(kN)': nom_load, 'Pu(kN)': ult_load, 'Mn(kN.m)': nom_moment, 'Mu(kN.m)': ult_moment, 'phi': phi_factor}
    df = pd.DataFrame(dict)

    sns.set_style('darkgrid')
    fig, ax = plt.subplots(figsize=(13,7))
    ax = sns.scatterplot(x='Mn(kN.m)', y='Pn(kN)', data=df, color='g')
    sns.scatterplot(x='Mu(kN.m)', y='Pu(kN)', data=df, color='r')
    ax.set_xlabel('Moment Capacity (kN.m)')
    ax.set_ylabel('Axial Load Capacity (kN)')
    plt.title('P-M Interaction Diagram')
    plt.show()

    phi = 0.65
    alpha = 0.85
    Ag = b*h
    As = math.pi * dbar**2 * ybar
    Pmax = phi*alpha*0.85*fc*(Ag-As) + As*fy /1000

    print('Maximum allowable concentric axial load: {} kN'.format(round(Pmax)))




class ColumnLoad:
    def __init__(self, ID, Label, Story, Load_Name, Station, P, V2, V3, M2, M3, T):
        self.ID = float(ID)
        self.Label = Label
        self.Story = Story
        self.Load_Name = Load_Name
        self.Station = float(Station)
        self.P = float(P)
        self.V2 = float(V2)
        self.V3 = float(V3)
        self.M2 = float(M2)
        self.M3 = float(M3)
        self.T = float(T)

    def __str__(self):
        return f"ColumnLoad(ID={self.ID}, Label='{self.Label}', Story='{self.Story}', Load_Name='{self.Load_Name}', Station={self.Station}, P={self.P}, V2={self.V2}, V3={self.V3}, M2={self.M2}, M3={self.M3}, T={self.T})"

class Column:
    def __init__(self, col_sect, col_type, fc, b, h, dbar, xbar, ybar, l, ke2, ke3, cover, fy, tie_db, tie_spc):
        self.col_sect = col_sect
        self.col_type = col_type
        self.fc = float(fc)
        self.b = float(b)
        self.h = float(h)
        self.dbar = float(dbar)
        self.xbar = float(xbar)
        self.ybar = float(ybar)
        self.l = float(l)
        self.ke2 = float(ke2)
        self.ke3 = float(ke3)
        self.cover = float(cover)
        self.tie_db = float(tie_db)
        self.tie_spc = float(tie_spc)

    def __str__(self):
        return f"Column(col_sect={self.col_sect}, col_type={self.col_type}, fc={self.fc}, b={self.b}, h={self.h}, dbar={self.dbar}, xbar={self.xbar}, ybar={self.ybar}, l={self.l}, ke2={self.ke2}, ke3={self.ke3}, cover={self.cover}, fy={self.fy}, tie_db={self.tie_db}, tie_spc={self.tie_spc})"


    from etabs_api import get_selected_columns

    # Retrieve selected columns from the ETABS model
    def selected_columns_data = get_selected_columns()

        # Create instances of the local Column class for each selected column
        columns = []
        for col_data in selected_columns_data:
            col_instance = Column(
                col_sect=col_data['section'],
                col_type=col_data['type'],
                fc=col_data['fc'],
                b=col_data['b'],
                h=col_data['h'],
                dbar=col_data['dbar'],
                xbar=col_data['xbar'],
                ybar=col_data['ybar'],
                l=col_data['length'],
                ke2=col_data['ke2'],
                ke3=col_data['ke3'],
                cover=col_data['cover'],
                fy=col_data['fy'],
                tie_db=col_data['tie_db'],
                    tie_spc=col_data['tie_spc']
            )
            columns.append(col_instance)
            print(f"Column added: {col_instance}")



C01 = Column('RECT', 'BRACED', 65, 900, 1100, 32, 4, 6, 3100, 1.0, 1.0, 40, 500, 0, 0)

#calculate_column_design("RECT", "BRACED", 65, 900, 1100, 32, 4, 6, 3100, 1.0, 1.0, 40, 500, 0, 0)