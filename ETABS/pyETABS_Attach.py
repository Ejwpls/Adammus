import os
import sys
import comtypes.client

# set the following flag to True to attach to an existing instance of the program
# otherwise a new instance of the program will be started
AttachToInstance = True

# set the following flag to True to manually specify the path to ETABS.exe
# this allows for a connection to a version of ETABS other than the latest installation
# otherwise the latest installed version of ETABS will be launched
SpecifyPath = False

# if the above flag is set to True, specify the path to ETABS below
ProgramPath = "C:\Program Files\Computers and Structures\ETABS 19\ETABS.exe"

if AttachToInstance:
    # attach to a running instance of ETABS
    try:
        # get the active ETABS object
        myETABSObject = comtypes.client.GetActiveObject("CSI.ETABS.API.ETABSObject")
    except (OSError, comtypes.COMError):
        print("No running instance of the program found or failed to attach.")
        sys.exit(-1)

else:
    # create API helper object
    helper = comtypes.client.CreateObject('ETABSv1.Helper')
    helper = helper.QueryInterface(comtypes.gen.ETABSv1.cHelper)
    if SpecifyPath:
        try:
            # 'create an instance of the ETABS object from the specified path
            myETABSObject = helper.CreateObject(ProgramPath)
        except (OSError, comtypes.COMError):
            print("Cannot start a new instance of the program from " + ProgramPath)
            sys.exit(-1)
    else:

        try:
            # create an instance of the ETABS object from the latest installed ETABS
            myETABSObject = helper.CreateObjectProgID("CSI.ETABS.API.ETABSObject")
        except (OSError, comtypes.COMError):
            print("Cannot start a new instance of the program.")
            sys.exit(-1)

# create SapModel object
SapModel = myETABSObject.SapModel
ret = SapModel.SetPresentUnits(6)

# get ModelFileName
ModelName = SapModel.GetModelFilename(False)

print(ModelName)
