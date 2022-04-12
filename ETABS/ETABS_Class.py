class frame:

    # Class Attribute
    # eltype = {"COLUMN","BEAM", "PIER", "SPANDREL")
    eltype = "COLUMN"

    # instance attribute
    def __init__(self,attributes,geometry,loads,reinforcement):
        self.attributes=[]