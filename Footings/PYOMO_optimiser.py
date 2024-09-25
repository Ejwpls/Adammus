import os
os.environ['PATH'] = r"C:\coin_solvers;" + os.environ['PATH']
import pyomo.environ as pyo

def optimize_foundation_pyomo(ColL, ColW, Pdl, Pll, BP):
    model = pyo.ConcreteModel()

    # Define the sets of allowed values (slightly expanded)
    L_values = [i/10 for i in range(60, 81, 1)]  # 6.0 to 8.0 in steps of 0.1
    D_values = [i/10 for i in range(11, 21, 1)]  # 1.1 to 2.0 in steps of 0.1
    fc_values = [25, 32, 40]  # Concrete strength options

    # Define index sets
    model.L_index = pyo.Set(initialize=range(len(L_values)))
    model.D_index = pyo.Set(initialize=range(len(D_values)))
    model.fc_index = pyo.Set(initialize=range(len(fc_values)))

    # Define binary variables
    model.L_binary = pyo.Var(model.L_index, domain=pyo.Binary)
    model.D_binary = pyo.Var(model.D_index, domain=pyo.Binary)
    model.fc_binary = pyo.Var(model.fc_index, domain=pyo.Binary)

    # Ensure only one value is chosen for each variable
    model.L_choice = pyo.Constraint(expr=sum(model.L_binary[i] for i in model.L_index) == 1)
    model.D_choice = pyo.Constraint(expr=sum(model.D_binary[i] for i in model.D_index) == 1)
    model.fc_choice = pyo.Constraint(expr=sum(model.fc_binary[i] for i in model.fc_index) == 1)

    # Define expressions for actual values
    model.L = sum(L_values[i] * model.L_binary[i] for i in model.L_index)
    model.W = model.L + (ColW - ColL)
    model.D = sum(D_values[i] * model.D_binary[i] for i in model.D_index)
    model.fc = sum(fc_values[i] * model.fc_binary[i] for i in model.fc_index)

    # Constraints
    # Bearing pressure constraint
    model.bearing_pressure = pyo.Constraint(expr=(Pdl + Pll + 6 * model.L * model.W * model.D) / (model.L * model.W) <= 1.00 * BP)

    # Objective function (minimize volume)
    model.cost = pyo.Objective(expr=model.L * model.W * (3-model.D)*model.D*(model.fc / 25 ), sense=pyo.minimize)

    # Solve
    solver = pyo.SolverFactory('couenne')
    solver.options['max_iter'] = 20000
    results = solver.solve(model, tee=True)

    if results.solver.status == pyo.SolverStatus.ok:
        return {
            'L': pyo.value(model.L),
            'W': pyo.value(model.W),
            'D': pyo.value(model.D),
            'cost': pyo.value(model.cost)
        }
    else:
        raise ValueError("Optimization failed.")

# Example usage
result = optimize_foundation_pyomo(ColL=0.8, ColW=0.5, Pdl=8000, Pll=1500, BP=250)
print("Optimal foundation parameters:")
for key, value in result.items():
    print(f"{key}: {value:.6f}")
