import pyETABS_Attach
import pandas as pd


def getDrift(loadcase):
    NumberResults = 0
    Story = []
    StepType = []
    StepNum = []
    Direction = []
    Drift = []
    Label = []
    LoadCase = []
    X = []
    Y = []
    Z = []

    ret = pyETABS_Attach.SapModel.Results.Setup.SetCaseSelectedForOutput(loadcase)
    ret = pyETABS_Attach.SapModel.Results.StoryDrifts(NumberResults, Story, LoadCase, StepType, StepNum, Direction,
                                                      Drift, Label, X, Y, Z)

    ret = pd.DataFrame(ret[1:-1]).T
    ret.columns = ['Story', 'LoadCase', 'StepType', 'StepNum', 'Direction', 'Drift', 'Label', 'X', 'Y', 'Z']
    # print(ret)

    return ret


print(getDrift("SPEC XY-DRIFT").head())
