from gekko import GEKKO
 
# Create a GEKKO model
m = GEKKO(remote=False)
 
# Define discrete sets for each variable
l_values = [1, 2, 3, 4, 5]  # Example values for length
w_values = [1, 2, 3, 4, 5]  # Example values for width
d_values = [1, 2, 3]        # Example values for depth
fc_values = [25, 32, 40, 50, 64]  # Example values for material cost factor
 
# Define variables using SOS (Special Ordered Sets)
l = m.sos1(l_values)
w = m.sos1(w_values)
d = m.sos1(d_values)
fc = m.sos1(fc_values)
 
# Define the objective function (e.g., minimize cost)
m.Minimize(l * w * d * fc)
 
# Add any necessary constraints (example constraint)
m.Equation(l + w + d <= 10)
 
# Set solver options to use APOPT
m.options.SOLVER = 1
 
# Solve the optimization problem
m.solve(disp=True)
 
# Display results
print(f"L: {l.VALUE[0]}")
print(f"W: {w.VALUE[0]}")
print(f"D: {d.VALUE[0]}")
print(f"fc: {fc.VALUE[0]}")