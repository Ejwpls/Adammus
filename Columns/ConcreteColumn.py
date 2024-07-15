import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


class ConcColumn:
    def __init__(self, sect, type="BRACED", fc=40, b=None, d=None, dbar=None, xbar=None, ybar=None, l=None, ke2=1, ke3=1, cover=30, tie_db=None, tie_spc=None):
        self.sect = sect
        self.type = type if type is not None else "BRACED"
        self.fc = float(fc) if fc is not None else 40 #[MPa] Concrete Strength
        self.b = float(b) #[mm] Column Width
        self.d = float(d) #[mm] Column Depth
        self.dbar = float(dbar) if dbar is not None else None #[mm] Rebar Diameter
        self.xbar = float(xbar) if xbar is not None else None #[No. of Bars]
        self.ybar = float(ybar) if ybar is not None else None #[No. of Bars]
        self.l = float(l) #[mm]
        self.ke2 = float(ke2) if ke2 is not None else 1
        self.ke3 = float(ke3) if ke3 is not None else 1
        self.cover = float(cover) if cover is not None else 30 #[mm]
        self.tie_db = float(tie_db) if tie_db is not None else None #[mm]
        self.tie_spc = float(tie_spc) if tie_spc is not None else None #[mm]
        
        #Constant Parameters
        self.fsy = 500                  # [MPa] Steel Yield Strength
        self.Es = 200000                # [MPa] Steel Modulus (Constant)
        self.ec = 0.003                 # [MPa] Concrete Maximum strain (MPa)
        self.ey = self.fsy / self.Es    # Steel strain at yield point
 

        es_squash = 0.0025
    
    def __str__(self):
        return f"ConcColumn(col_sect='{self.sect}', col_type='{self.type}', fc={self.fc}, b={self.b}, d={self.d}, dbar={self.dbar}, xbar={self.xbar}, ybar={self.ybar}, l={self.l}, ke2={self.ke2}, ke3={self.ke3}, cover={self.cover}, tie_db={self.tie_db}, tie_spc={self.tie_spc})"
    
    #Column Geometric Properties
    @property
    def Agross(self):
        if self.sect == 'RECT':
            return self.b * self.d
        elif self.sect == 'CIRC':
            return math.pi * (self.d/2) ** 2
        else:
            raise ValueError("Invalid section type")
        
    # Moment of Inertia [mm4]
    def Ig(self, sect, b, d):
        if sect == 'RECT':
            return b * d ** 3 / 12
        elif sect == 'CIRC':
            return math.pi / 4 * (d/2) ** 4
    
    @property
    def Ig3(self):
        return self.Ig(self.sect, self.b, self.d)
    
    @property
    def Ig2(self):
        return self.Ig(self.sect, self.d, self.b)
    
    #Radius of Gyration [mm]
    def rgy(self,Ig,Agross):
        return (Ig / Agross) ** 0.5

    @property
    def rgy3(self):
        return self.rgy(self.Ig3,self.Agross)
    
    @property
    def rgy2(self):
        return self.rgy(self.Ig2,self.Agross)
    
    #Reinforcement Arrangement 
    @property
    def barspc3(self):
        return self.d03 / (self.ybar -1)
    
    @property
    def barspc2(self):
        return self.d02 / (self.xbar -1)
        
    def bar_arrangement(self):
        # Column Rebar - Rect Section
        i_y = []
        i_yz = []
        i_x = []
        i_xz = []

        xbar = int(self.xbar)
        ybar = int(self.ybar)

        # No. Bars in y_direction
        i_y = [xbar]
        i_yz = [self.d03/2]
        for i in range(1, ybar-1):
            # Vertical layer
            i_y.append(2)
            i_yz.append(self.d03/2 - self.barspc3*i)
        i_y.append(xbar)
        i_yz.append(-self.d03/2)

        # No. Bars in x_direction
        i_x = [ybar]
        i_xz = [self.d02/2]
        for i in range(1, xbar-1):
            # Horizontal layer
            i_x.append(2)
            i_xz.append(self.d02/2 - self.barspc2*i)
        i_x.append(ybar)
        i_xz.append(-self.d02/2)

        # Assign coordinate for each bar
        xi = []
        yi = []
        marker_size = []

        for k in range(0, ybar):
            for j in range(0, i_y[k]):
                if i_y[k] == 2 and j == 1:
                    xi.append(i_xz[xbar-1])
                    yi.append(i_yz[k])
                else:
                    xi.append(i_xz[j])
                    yi.append(i_yz[k])

                marker_size.append(self.dbar)

        return xi, yi, marker_size
    
    #function to plot bar arrangemnet in the column
    def plot_bar_arrangement(self):
        xi, yi, marker_size = self.bar_arrangement()
        scaled_marker_size = np.array(marker_size) * 5  # Scale the marker size by a factor of 10

        plt.figure(figsize=(8, 8))  # Set the figure size to 8x8 inches
        plt.scatter(xi, yi, s=scaled_marker_size)
        plt.xlabel('X-axis (mm)')
        plt.ylabel('Y-axis (mm)')
        plt.title('Column Bar Arrangement')

        # Calculate the center of the column
        center_x = self.b / 2
        center_y = self.d / 2

        # Calculate the buffer size
        buffer = 100

        # Set the axis limits with buffer
        plt.xlim(-center_x - buffer, center_x + buffer)
        plt.ylim(-center_y - buffer, center_y + buffer)

        # Plot column outline with the center as the origin
        plt.plot([-center_x, center_x, center_x, -center_x, -center_x], [-center_y, -center_y, center_y, center_y, -center_y], 'k-', linewidth=2)

        plt.show()

    
    #Distances
    def deff(self,dirn):
        if dirn == '3':
            return self.d - self.cover - self.tie_db - self.dbar/2
        elif dirn == '2':
            return self.b - self.cover - self.tie_db - self.dbar/2
    
    @property
    def d03(self):
        return self.d - self.cover - self.tie_db - self.dbar/2
    
    @property
    def d02(self):
        return self.b - self.cover - self.tie_db - self.dbar/2
    
    @property
    def dprime(self):
        return self.cover + self.tie_db + self.dbar/2
    
    @property
    def total_bars(self):
        if self.sect == 'CIRC':
            return self.ybar
        elif self.sect == 'RECT':
            return 2*(self.xbar + self.ybar - 2)
    
    def total_Ast(self):
        return math.pi * self.dbar**2 / 4 * self.total_bars()
    
    #Concrete Parameters - CL 10.6
    @property
    def alpha1(self):
        return min(0.85, max(1 - 0.003 * self.fc, 0.72))       # [MPa] Concrete Strength Reduction Factor

    @property
    def alpha2(self):
        alpha2 = max(0.85 - 0.0015 * self.fc, 0.67)     # [MPa] Concrete Strength Reduction Factor
        if self.sect == 'CIRC':
            alpha2 = alpha2 * 0.95
        return alpha2
    
    @property
    def gamma(self):
        return max(0.97-0.0025 * self.fc, 0.67)         # [MPa] Concrete Strength Reduction Factor
    
    #Squash Load
    def Nuo(self):
        Anet = (self.Agross()) - (self.total_Ast())
        Nuo = self.alpha1() * self.fc * Anet /1000 + (self.total_Ast() * self.fsy /1000)
        return Nuo
    
    #Maximum Axial tension Load
    def Ntens(self):
        return self.total_Ast() * self.fsy /1000
    

    #Decompression Point
    def Ndecomp(self):
        return self.alpha2() * self.fc * self.Agross() /1000
    

    def calculate_c(self, width, step):
        x= step
        c_list = []
        while x < (width-step):
            c_list.append(x)
            x += step
        return c_list
    
    def calculate_cb(self, d):
        return d*self.materials.concrete.Ɛu/(self.materials.concrete.Ɛu + self.materials.steel.Ɛy)
    
    def calculate_a(self, c, beta1):
        return c*beta1
    
    def calculate_d(self, width, coordinate):
        return width-coordinate

    def calculate_strain(self, c,d):
        return self.materials.concrete.Ɛu*(d-c)/c
    
    def calculate_stress(self, strain):
        stress = self.materials.steel.elasticity_modulus*strain
        if stress > self.materials.steel.fy:
            return self.materials.steel.fy
        elif stress < -self.materials.steel.fy:
            return -self.materials.steel.fy
        else:
            return stress
    
    def calculate_steel_force(self, stress, area):
        return stress * area/1000
    
    def calculate_steel_momentum(self, force, d,width):
        if self.geometry.units == "IMP":
            return force *(d-width/2)/12
        else:
            return force *(d-width/2)/1000
    
    
    def calculate_FC(self, a,base):
        return 0.85*self.materials.concrete.fc*a*base/1000
    
    def calculate_FC_momentum(self, a, base, width):
        if self.geometry.units == "IMP":
            return self.calculate_FC(a, base)*(width/2-a/2)/12
        else:
            return self.calculate_FC(a, base)*(width/2-a/2)/1000
        
    def calculate_ACI_phi(self, strain):
        if strain < self.materials.steel.Ɛy:
            return 0.65
        elif strain > self.materials.steel.Ɛy +0.003:
            return 0.90
        else:
            return 0.65+0.25*(strain-self.materials.steel.Ɛy)/(0.003)  

    def calculate_PM_values(self, c, width, coordinates, base, areas):
        d = [self.calculate_d(width, coordinate) for coordinate in coordinates]
        a = self.calculate_a(c, self.materials.concrete.beta_1)
        FC = self.calculate_FC(a, base)
        strains = [self.calculate_strain(c,x) for x in d]
        stresses = [self.calculate_stress(strain) for strain in strains]
        steel_forces = [self.calculate_steel_force(stress,area) for (stress,area) in zip(stresses, areas)]
        sum_steel_forces = sum(steel_forces)
        Pn= (FC - sum_steel_forces)
        if Pn > self.Pn_max: self.Pn_max
        if Pn < self.Pn_t: self.Pn_t
        FC_momentum = self.calculate_FC_momentum(a,base, width)
        steel_momentums = ([self.calculate_steel_momentum(force,x,width) for (force,x) in zip(steel_forces, d)])
        sum_steel_momentums = sum(steel_momentums)
        Mn = (FC_momentum + sum_steel_momentums)
        return (Mn, Pn)
    
    def calculate_PM_list(self, width, coordinates, base, areas, step):
        c_list = self.calculate_c(width, step)
        d = max(coordinates)
        cb = self.calculate_cb(d)
        c_list.append(cb)
        c_list.sort()
        PM_list = [self.calculate_PM_values(c,width,coordinates, base,areas) for c in c_list]
        PM_list.insert(0,(0,self.Pn_t))
        PM_list.append((0,self.Pn_max))

        return PM_list


    def calculate_PM_dataframe(self, width, coordinates, base, areas, step):
        d = [self.calculate_d(width, coordinate) for coordinate in coordinates]
        p_m_values = pd.DataFrame()
        c_list = self.calculate_c(width, step)
        cb = self.calculate_cb(max(d))
        c_list.append(cb)
        c_list.sort()
        a_list = [self.calculate_a(value, self.materials.concrete.beta_1) for value in c_list]
        FC_list = [self.calculate_FC(value, base) for value in a_list]
        p_m_values.insert(0, "c", c_list)
        p_m_values.insert(1,"a", a_list)
        p_m_values.insert(2, "FC", FC_list)
        d_columns = {f'd{i+1}': value for i, value in enumerate(d)}
        p_m_values = p_m_values.assign(**d_columns)
        #strains
        for i, value in enumerate(d):
            p_m_values[f'strain{i+1}'] = p_m_values["c"].apply(lambda x: self.calculate_strain(x,value))
            p_m_values[f'stress{i+1}'] = p_m_values[f'strain{i+1}'].apply(lambda x: self.calculate_stress(x))
            p_m_values[f'steel_force{i+1}'] = p_m_values[f'stress{i+1}'].apply(lambda x: self.calculate_steel_force(x,areas[i]))
        p_m_values['SUM_steel_force'] = p_m_values.filter(like="steel_force").sum(axis=1)
        p_m_values['Pn'] = np.where(p_m_values['FC'] - p_m_values['SUM_steel_force'] > self.Pn_max, self.Pn_max, np.where(p_m_values['FC'] - p_m_values['SUM_steel_force'] < self.Pn_t, self.Pn_t, p_m_values['FC'] - p_m_values['SUM_steel_force']))
        p_m_values['FC_momentum'] = p_m_values['a'].apply(lambda x: self.calculate_FC_momentum(x,base, width))
        for i, value in enumerate(d):
            p_m_values[f'steel_momentum{i+1}'] = p_m_values[f'steel_force{i+1}'].apply(lambda x: self.calculate_steel_momentum(x,value,width))
        p_m_values['SUM_steel_momentum'] = p_m_values.filter(like="steel_momentum").sum(axis=1)
        p_m_values['Mn'] = p_m_values['SUM_steel_momentum'] + p_m_values['FC_momentum']
        p_m_values['max_strain'] = p_m_values.filter(like="strain").max(axis=1)
        p_m_values['phi'] = p_m_values['strain1'].apply(lambda x: self.calculate_ACI_phi(x))
        p_m_values['phiPn'] = p_m_values['Pn']*p_m_values['phi']
        p_m_values['phiMn'] = p_m_values['Mn']*p_m_values['phi']
        return p_m_values
    
    @property
    def PM_table_1(self):
        return self.calculate_PM_dataframe(self.geometry.width1, self.positions_1, self.geometry.width2, self.areas_d1, 25)

    @property
    def PM_table_2(self):
        return self.calculate_PM_dataframe(self.geometry.width2, self.positions_2, self.geometry.width1, self.areas_d2, 25)

    @property
    def PM_direction_1_df(self):
        values = self.calculate_PM_dataframe(self.geometry.width1, self.positions_1, self.geometry.width2, self.areas_d1, 0.1)
        PM_list = list(zip(values['Mn'], values['Pn']))
        PM_list.insert(0,(0,self.Pn_t))
        PM_list.append((0,self.Pn_max))
        return PM_list
    
    @property
    def phi_PM_direction_1_df(self):
        values = self.calculate_PM_dataframe(self.geometry.width1, self.positions_1, self.geometry.width2, self.areas_d1, 0.1)
        PM_list = list(zip(values['phiMn'], values['phiPn']))
        PM_list.insert(0,(0,self.Pn_t*0.9))
        PM_list.append((0,self.Pn_max*0.65))
        return PM_list
    
    @property
    def PM_direction_1(self):
        return self.calculate_PM_list(self.geometry.width1, self.positions_1, self.geometry.width2, self.areas_d1, 0.1)
    
    @property
    def PM_direction_2_df(self):
        values = self.calculate_PM_dataframe(self.geometry.width2, self.positions_2, self.geometry.width1, self.areas_d2, 0.1)
        PM_list = list(zip(values['Mn'], values['Pn']))
        PM_list.insert(0,(0,self.Pn_t))
        PM_list.append((0,self.Pn_max))
        return PM_list
    
    @property
    def phi_PM_direction_2_df(self):
        values = self.calculate_PM_dataframe(self.geometry.width2, self.positions_2, self.geometry.width1, self.areas_d2, 0.1)
        PM_list = list(zip(values['phiMn'], values['phiPn']))
        PM_list.insert(0,(0,self.Pn_t*0.9))
        PM_list.append((0,self.Pn_max*0.65))
        return PM_list
    @property
    def PM_direction_2(self):
        return self.calculate_PM_list(self.geometry.width2, self.positions_2, self.geometry.width1, self.areas_d2, 0.1)
    
    def plot_PM(self, points1 = None, points2 = None):
        x_values1, y_values1 = zip(*self.PM_direction_1_df)
        # x_values11, y_values11 = zip(*self.phi_PM_direction_1_df)


        x_values2, y_values2 = zip(*self.PM_direction_2_df)
        # x_values22, y_values22 = zip(*self.phi_PM_direction_2_df)


        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(8, 6))
        

        ax1.plot(x_values1, y_values1, marker='o', linestyle='-', color = "#5AC990", label = "Direction 1")
        # ax1.plot(x_values11, y_values11, marker='o', linestyle='-', color = 'red', label = "phi Direction 1")
        ax2.plot(x_values2, y_values2, marker='x', linestyle='-', color = "#5AC990", label = "Direction 2")
        # ax2.plot(x_values22, y_values22, marker='x', linestyle='-', color = 'red', label = "phi Direction 2")

        ##Added points
        if points1:
            x_values_points, y_values_points = zip(*points1)
            # print(x_values_points, y_values_points)
            ax1.scatter(x_values_points, y_values_points, marker = 'x', color='blue', label='Pu, Mu1')
        
        if points2:
            x_values_points, y_values_points = zip(*points2)
            ax2.scatter(x_values_points, y_values_points, marker = 'o' ,color='blue', label='Pu, Mu2')
        #ax1.title("P-M Diagram Direction 1")
        if self.geometry.units == "IMP":
            ax1.set_xlabel('Mn1 [kips-ft]')
            ax1.set_ylabel('Pn1 [kips]')
        else:
            ax1.set_xlabel('Mn1 [kN-m]')
            ax1.set_ylabel('Pn1 [kN]')
        #ax1.grid(True, linestyle='--', alpha=0.5)       
        ax1.spines['top'].set_linestyle('-')
        ax1.spines['bottom'].set_linestyle('-')
        ax1.spines['left'].set_linestyle('-')
        ax1.spines['right'].set_linestyle('-')
        ax1.grid(True, linestyle='--', alpha=0.5)  # Set grid for the first subplot
        #ax1.legend("Direction 1")
 

        #ax2.title("P-M Diagram Direction 2")
        ax2.set_xlabel('Mn2 [kips-ft]')
        ax2.set_ylabel('Pn2 [kips]')
        #ax2.grid(True, linestyle='--', alpha=0.5)       
        ax2.spines['top'].set_linestyle('-')
        ax2.spines['bottom'].set_linestyle('-')
        ax2.spines['left'].set_linestyle('-')
        ax2.spines['right'].set_linestyle('-')
        ax2.grid(True, linestyle='--', alpha=0.5)  # Set grid for the first subplot
        #ax2.legend("Direction 2")
       
        # Create a single legend for both subplots and place it in the middle
        lines, labels = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        lines.extend(lines2)
        labels.extend(labels2)
        fig.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, 0.95))

        plt.show()
    # Plotting the DataFrame as a table using Matplotlib
    def plot_PM_table(self, table):
        fig, ax = plt.subplots(figsize=(6, 2))  # Adjust the figsize as needed
        ax.axis('tight')
        ax.axis('off')
        ax.table(cellText=table.values, colLabels=table.columns, loc='center')

        plt.show()


    def plot_PM_table_1(self):
        return self.plot_PM_table(self.PM_table_1[['Mn','Pn']])
    
    def plot_PM_table_2(self):
        return self.plot_PM_table(self.PM_table_2[['Mn','Pn']])
    
    



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


    #from etabs_api import get_selected_columns

    # Retrieve selected columns from the ETABS model
    #selected_columns_data = get_selected_columns()

    # Create instances of the local Column class for each selected column
    