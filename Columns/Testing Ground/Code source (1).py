from math import pi, sqrt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches


class Geometry:
    def __init__(self, width1: float = 10.0, width2: float = 12, width3: float = 15.0, units: str = "IMP"):
        self.units = units
        self.width1 = width1
        self.width2 = width2
        if self.units == "IMP":
            self.width3 = width3 * 12
        else:
            self.width3 = width3 * 1000

    #Geometry Properties    
    @property
    def cross_section_area12(self):
        return self.width1 * self.width2

    @property
    def cross_section_area13(self):
        return self.width1 * self.width3
    
    @property
    def cross_section_area32(self):
        return self.width3 * self.width2
    
    #Geometry Coordinates
    @property
    def coordinates(self):
        return [(0,0),(0, self.width2), (self.width1, self.width2), (self.width1, 0),(0,0)]


class Rebar:
    def __init__(self, As: int = 5, units: str = "IMP"):
        self.As = As
        self.units = units

    @property
    def bar_diameter(self) -> float:
        """Get the diameter of the bar in inches."""
        if self.units == "IMP":
            return self.As / 8
        else:
            return self.As
    
    @property
    def bar_area(self) -> float:
        """Get the cross-sectional area of the bar in square inches."""
        return (self.bar_diameter ** 2) * pi / 4


class Longitudinal_reinforcement:

    def __init__(self, n1: int = 3, As1: Rebar = Rebar(5), n2: int = 4, As2: Rebar = Rebar(5), As_corner: Rebar = Rebar(5)):
        self.n1 = n1
        self.As1 = As1
        self.n2 = n2
        self.As2 = As2
        self.As_corner = As_corner
    
    @property
    def Ast(self):
        return (self.As_corner.bar_area * 4 + self.As1.bar_area * (self.n1 - 2) * 2 + self.As2.bar_area * (self.n2 - 2) * 2)
    
    @property
    def nAs(self):
        return self.n1*2 + (self.n2-2)*2


class Transverse_reinforcement:
    def __init__(self, legs1: int = 3, legs2: int = 3, As_t: Rebar = Rebar(3), s1: float = 8, s2: float = 10):
        self.legs1 = legs1
        self.legs2 = legs2
        self.As_t = As_t
        self.s1 = s1
        self.s2 = s2
    
    @property
    def Ast_1(self):
        return self.legs1*self.As_t.bar_area
    
    @property
    def Ast_11(self):
        return self.Ast_1/self.s1
    
    @property
    def Ast_12(self):
        return self.Ast_1/self.s2
    
    @property
    def Ast_2(self):
        return self.legs2*self.As_t.bar_area
    
    @property
    def Ast_21(self):
        return self.Ast_2/self.s1
    
    @property
    def Ast_12(self):
        return self.Ast_2/self.s2


class Concrete:

    def __init__(self, fc: int = 4000, cover: float = 1.5, Ɛu: float = 0.003, units: str = "IMP", EC_factor: float = 57000) -> None:
        self.fc = fc
        self.cover = cover
        self.Ɛu = Ɛu
        self.units = units
        self.EC_factor = EC_factor
    
    @property
    def elastic_modulus(self) -> float:
        return self.EC_factor * sqrt(self.fc)

    @property
    def beta_1(self):
        if self.units == "IMP":
            conversion_factor = 1
        else:
            conversion_factor = .007

        if self.fc <= 4000 * conversion_factor:
            return 0.85
        elif self.fc <= 5000 * conversion_factor:
            return 0.80
        elif self.fc <= 6000 * conversion_factor:
            return 0.75
        elif self.fc <= 7000 * conversion_factor:
            return 0.70
        else:
            return 0.65


class Steel:
    def __init__(self, fy: float = 60000, Es: float = 29000000):
        self.fy = fy
        self.Es = Es
    
    @property
    def elasticity_modulus(self) -> float:
        return self.Es

    @property
    def Ɛy(self):
        return self.fy / self.elasticity_modulus


class Materials:
    def __init__(self, concrete: Concrete = Concrete(), steel: Steel = Steel()):
        self.concrete = concrete
        self.steel = steel
    
    @property
    def modular_factor(self) -> float:
        return self.steel.elasticity_modulus / self.concrete.elastic_modulus

    
class CR_Columns:
    def __init__(self, geometry: Geometry = Geometry(), longitudinal: Longitudinal_reinforcement  = Longitudinal_reinforcement(), 
                 transverse: Transverse_reinforcement = Transverse_reinforcement(), materials: Materials = Materials(), code = None):
        self.geometry = geometry
        self.longitudinal = longitudinal
        self.transverse = transverse
        self.materials = materials
        self.first_position = self.materials.concrete.cover + self.transverse.As_t.bar_diameter + self.longitudinal.As_corner.bar_diameter/2
        self.last_position1 = self.geometry.width1-self.first_position
        self.last_position2 = self.geometry.width2-self.first_position
    
    #Initial methods
    def set_space(self, width):
        return width - self.first_position*2
    
    def set_spacing(self, width, n):
        return width/(n-1)
    
    def set_free_spacing(self, diameter, spacing):
        return spacing + diameter
    
    def set_positions(self, n, spacing):
        x = self.first_position
        positions = []
        counter = 0
        while counter < n:
            positions.append(x)
            x += spacing
            counter += 1
        return positions

    #Properties
    @property
    def free_width1(self):
        return self.set_space(self.geometry.width1)

    @property
    def free_width2(self):
        return self.set_space(self.geometry.width2)
    
    #Viability
    @property
    def viability1(self):
        return self.free_width1 > self.transverse.As_t.bar_diameter*6
    
    @property
    def viability2(self):
        return self.free_width2 > self.transverse.As_t.bar_diameter*6
    
    @property
    def viability(self):
        return self.viability1 & self.viability2
    
    #Reinforcement Properties

    @property
    def stirrups_coordinates(self):
        initial = self.materials.concrete.cover + self.transverse.As_t.bar_diameter/2
        return [
            (initial,initial),(initial, self.geometry.width2-initial), 
            (self.geometry.width1-initial, self.geometry.width2-initial), 
            (self.geometry.width1-initial, initial),(initial,initial)
            ]
    @property
    def spacing_1(self):
        return self.set_spacing(self.free_width1,self.longitudinal.n1)
    
    @property
    def spacing_2(self):
        return self.set_spacing(self.free_width2,self.longitudinal.n2)
    
    @property
    def free_spacing_1(self):
        return self.set_free_spacing(self.longitudinal.As1.bar_diameter, self.spacing_1)
    
    @property
    def free_spacing_2(self):
        return self.set_free_spacing(self.longitudinal.As2.bar_diameter, self.spacing_2)
    
    @property
    def positions_1(self):
        return self.set_positions(self.longitudinal.n1, self.spacing_1)

    @property
    def positions_2(self):
        return self.set_positions(self.longitudinal.n2, self.spacing_2)

    
    @property
    def reinforcement_coordinates(self):
        reinforcement_list = {
            "Corner_Bars":{
                1:[
                (self.first_position,self.first_position),
                (self.last_position1, self.first_position)],
                2:[ 
                (self.last_position1, self.last_position2),
                (self.first_position, self.last_position2)]},
            "Direction_1":{
                1:[(x,self.first_position) for x in self.positions_1[1:-1]],
                2:[(x,self.last_position2) for x in self.positions_1[1:-1]]
                           },
            "Direction_2":{
                1:[(self.first_position,y) for y in self.positions_2[1:-1]],
                2:[(self.last_position1,y) for y in self.positions_2[1:-1]]
                           },
            "Direction_3":{1:[],},
        }
        return reinforcement_list
    
    @property
    def areas_d1(self):
        corners = self.longitudinal.As_corner.bar_area * 2 + self.longitudinal.As2.bar_area * (self.longitudinal.n2-2)
        areas = [corners, ]
        n = 0
        while n < self.longitudinal.n1-2:
            areas.append(self.longitudinal.As1.bar_area*2)
            n +=1
        areas.append(corners)
        return areas

    @property
    def areas_d2(self):
        corners = self.longitudinal.As_corner.bar_area * 2 + self.longitudinal.As1.bar_area * (self.longitudinal.n1-2)
        areas = [corners, ]
        n = 0
        while n < self.longitudinal.n2 - 2:
            areas.append(self.longitudinal.As2.bar_area*2)
            n += 1
        areas.append(corners)
        return areas 
    
    @property
    def stirrup_plot_coordinates(self):
        initial = self.materials.concrete.cover + self.transverse.As_t.bar_diameter*3.5 + self.transverse.As_t.bar_diameter/2
        initial2 = self.materials.concrete.cover + self.transverse.As_t.bar_diameter/2
        return [(initial2, initial), (initial2, self.geometry.width2-initial), 
                (initial, self.geometry.width2-initial2), (self.geometry.width1 - initial, self.geometry.width2-initial2),
                (self.geometry.width1-initial2, self.geometry.width2-initial), (self.geometry.width1-initial2, initial),
                (self.geometry.width1-initial, initial2), (initial,initial2)
                ]
    #Internal Methods
    ##Setter Methods
    ###Longitudinal Reinforcement
    def set_longitudinal_rebar(self, As: Rebar, coordinates, direction):
        bar_dict={
            "As": As,
            "Coordinates": tuple(coordinates),
        }
        self.reinforcement[direction].append()


    def plot_cross_section(self):
        if self.viability:
            fig, ax = plt.subplots()
            n=0
            diameters = [self.longitudinal.As_corner.bar_diameter, self.longitudinal.As1.bar_diameter, self.longitudinal.As2.bar_diameter, self.transverse.As_t.bar_diameter]
            ## Plot Rectangle
            for x, y in self.geometry.coordinates:
                ax.plot(x, y, 'grey')
                    
            rectangle = patches.Rectangle((0 ,0), self.geometry.width1,self.geometry.width2,
                                        fill=True, facecolor="white", edgecolor="#6B7280")       
            ax.add_patch(rectangle)
            #arrow(0,-2,self.free_width1/2,0)
            ## Plot Rebars
            for value in self.reinforcement_coordinates.values():
                # Plot the reinforcement bars
                for l in value.values():
                    for x, y in l:
                        # Create a circle patch for each reinforcement bar
                        circle = patches.Circle((x, y), diameters[n]/2, fill=True, facecolor= "#A3A8B0", edgecolor="#A3A8B0")
                        ax.add_patch(circle)
                n += 1
            ## Plot stirrups
            if self.geometry.units == "IMP":
                linewidth = self.transverse.As_t.bar_diameter*10
            else:
                linewidth = self.transverse.As_t.bar_diameter/2

            initial = self.materials.concrete.cover + self.transverse.As_t.bar_diameter*3.5 + self.transverse.As_t.bar_diameter/2
            archs_centers = [(initial, initial), (initial, self.geometry.width2-initial), 
                            (self.geometry.width1-initial, self.geometry.width2-initial),
                            (self.geometry.width1-initial, initial)
                            ]
            for i in range(0, len(self.stirrup_plot_coordinates), 2):
                x1, y1 = self.stirrup_plot_coordinates[i]
                x2, y2 = self.stirrup_plot_coordinates[i + 1]
                    
                # Plot the line
                ax.plot([x1, x2], [y1, y2], '#6B7280',linewidth=linewidth)

            diameter = self.transverse.As_t.bar_diameter * 7
                
            #Draw arcs in the corners.

            for i in range(len(archs_centers)):
                # Extract coordinates for the current iteration
                center_x, center_y = archs_centers[i]
                arc = patches.Arc((center_x, center_y), width=diameter, height=diameter, theta1=180-i*90, theta2=270-i*90, color='#6B7280', linewidth=linewidth)
                ax.add_patch(arc)

            # Set labels and title
            if self.geometry.units == "IMP":
                units = "[in]"
                
            else:
                units = "[mm]"

            ax.set_xlabel(f"{self.geometry.width1} {units}")
            ax.set_ylabel(f"{self.geometry.width2} {units}")        
            ax.set_title("Cross-Section")
            ax.set_aspect('equal')               
            ax.set_xticks([])
            ax.set_yticks([])
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            # Show the plot
            plt.show()
        else:
            # Create a figure and axis
            fig, ax = plt.subplots()

            # Define rectangle dimensions
            width = 400
            height = 500

            # Create a rectangle patch
            rectangle = patches.Rectangle((0, 0), width, height, linewidth=1, edgecolor='#6B7280', facecolor='none')

            # Add the rectangle to the plot
            ax.add_patch(rectangle)

            # Add text to the center of the rectangle
            text = "Figure not available"
            text_size = 14

            text_x = width / 2
            text_y = height / 2

            ax.text(text_x, text_y, text, fontsize=text_size, color='#6B7280', ha='center', va='center')

            # Set axis limits
            ax.set_xlim(0, width)
            ax.set_ylim(0, height)

            # Remove axes
            ax.axis('off')

            # Show the plot
            plt.show()

        
    #P-M Interaction Diagram
    @property
    def Pn_max(self):
        return (0.80*0.85*self.materials.concrete.fc*(self.geometry.cross_section_area12-self.longitudinal.Ast)+self.materials.steel.fy*self.longitudinal.Ast)/1000
    
    @property
    def Pn_t(self):
        return -self.longitudinal.Ast * self.materials.steel.fy/1000

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

width1 = _page_width_1.to(ct.units("in"))
width2 = _page_width_2.to(ct.units("in"))
fc = _page_f_c.to(ct.units("psi"))
cover = _page_clear_cover.to(ct.units("in"))
fy = _page_fy.to(ct.units("psi"))



geometry = Geometry(width1.magnitude, width2.magnitude)
longitudinal = Longitudinal_reinforcement(_page_N1_bars, Rebar(_page_bar1_diameter), _page_N2_bars, Rebar(_page_bar2_diameter),Rebar(_page_corner_bar_diameter))
transversal = Transverse_reinforcement(2,2,Rebar(_page_Stirrup_diameter),8,10)
concrete = Concrete(fc.magnitude, cover.magnitude)
steel = Steel(fy.magnitude)
materials = Materials(concrete, steel)
column = CR_Columns(geometry, longitudinal, transversal, materials)
pm = column.plot_PM()